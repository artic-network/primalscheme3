import hashlib
import json
import pathlib
import shutil
import sys

from primaldimer_py import do_pools_interact_py  # type: ignore

# Core imports
from primalscheme3.__init__ import __version__
from primalscheme3.core.bedfiles import BedPrimerPair, read_in_bedprimerpairs
from primalscheme3.core.digestion import (
    DIGESTION_ERROR,
    f_digest_to_count,
    r_digest_to_count,
)
from primalscheme3.core.logger import setup_loger
from primalscheme3.core.msa import MSA
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.core.thermo import THERMORESULT, forms_hairpin, thermo_check_kmers


class BedPrimerPairRepair(BedPrimerPair):
    def __init__(self, BedPrimerPair):
        super().__init__(
            fprimer=BedPrimerPair.fprimer,
            rprimer=BedPrimerPair.rprimer,
            msa_index=BedPrimerPair.msa_index,
            chrom_name=BedPrimerPair.chrom_name,
            amplicon_prefix=BedPrimerPair.amplicon_prefix,
            amplicon_number=BedPrimerPair.amplicon_number,
            pool=BedPrimerPair.pool,
        )
        self.new_fkmer_seqs = set()
        self.new_rkmer_seqs = set()

    def to_bed(self, new_primer_prefix: str) -> str:
        """
        Custom to bed function to handle the new seqs
        """
        num_fprimers = len(self.fprimer.seqs)
        num_rprimers = len(self.new_fkmer_seqs)

        current_fprimers_str = self.fprimer.__str__(
            reference=f"{self.chrom_name}",
            amplicon_prefix=f"{self.amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        )

        new_fprimer_str = [
            f"{self.chrom_name}\t{self.fprimer.end-len(x)}\t{self.fprimer.end}\t{new_primer_prefix}_{self.amplicon_number}_LEFT_{count}\t{self.pool + 1}\t+\t{x}"
            for count, x in enumerate(self.new_fkmer_seqs, start=num_fprimers)
        ]

        current_rpimer_str = self.rprimer.__str__(
            reference=f"{self.chrom_name}",
            amplicon_prefix=f"{self.amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        )

        new_rprimer_str = [
            f"{self.chrom_name}\t{self.rprimer.start}\t{self.rprimer.start + len(x)}\t{new_primer_prefix}_{self.amplicon_number}_RIGHT_{count}\t{self.pool + 1}\t+\t{x}"
            for count, x in enumerate(self.new_rkmer_seqs, start=num_rprimers)
        ]

        return "\n".join(
            [
                x.strip()
                for x in [
                    current_fprimers_str,
                    *new_fprimer_str,
                    current_rpimer_str,
                    *new_rprimer_str,
                ]
            ]
        )


class SeqStatus:
    seq: str | None
    count: int
    thermo_status: THERMORESULT | DIGESTION_ERROR

    def __init__(
        self, seq: str | None, count: int, thermo_status: THERMORESULT | DIGESTION_ERROR
    ):
        self.seq = seq
        self.count = count
        self.thermo_status = thermo_status


def detect_early_return(seq_counts: dict[str | DIGESTION_ERROR, int]) -> bool:
    """
    Checks for an early return condition, will return True condition is met
    """
    # Check for early return conditions
    for error, count in seq_counts.items():
        if count == -1 and type(error) == DIGESTION_ERROR:
            return True
    return False


def thermo_check_seq(seq, base_cfg: dict) -> THERMORESULT | DIGESTION_ERROR:
    # Thermo check
    themo_result = thermo_check_kmers([seq], base_cfg)
    if themo_result != THERMORESULT.PASS:
        return themo_result

    # Check for hairpins
    elif forms_hairpin([seq], base_cfg):  # type: ignore
        return DIGESTION_ERROR.HAIRPIN_FAIL

    # Check for Homodimer
    if do_pools_interact_py([seq], [seq], base_cfg["dimerscore"]):
        return DIGESTION_ERROR.DIMER_FAIL

    return THERMORESULT.PASS


def report_check(
    seqstatus: SeqStatus,
    current_primer_seqs: set[str],
    seqs_in_pools: list[list[str]],
    pool: int,
    dimerscore: float,
    logger,
) -> bool:
    """
    Will carry out the checks and report the results via the logger. Will return False if the seq should not be added
    """

    report_seq = seqstatus.seq if seqstatus.seq is not None else "DIGESTION_ERROR"
    report_seq = report_seq.rjust(40, " ")

    # Check it passed thermo
    if seqstatus.thermo_status != THERMORESULT.PASS or seqstatus.seq is None:
        logger.warning(
            "{newseq}\t{count}\t<red>FAILED</>: {thermo_status}",
            thermo_status=seqstatus.thermo_status,
            newseq=report_seq,
            count=seqstatus.count,
        )
        return False

    # Check it is a new seq
    if seqstatus.seq in current_primer_seqs:
        logger.warning(
            "{newseq}\t{count}\t<blue>NOTADDED</>: Already Included",
            newseq=report_seq,
            count=seqstatus.count,
        )
        return False

    # Check for minor allele
    if seqstatus.count < 0:
        logger.warning(
            "{newseq}\t{count}\t<red>FAILED</>: Minor Allele",
            newseq=report_seq,
            count=seqstatus.count,
        )
        return False

    # Check for dimer with pool
    if do_pools_interact_py([seqstatus.seq], seqs_in_pools[pool], dimerscore):
        logger.warning(
            "{newseq}\t{count}\t<red>FAILED</>: Interaction with pool",
            newseq=report_seq,
            count=seqstatus.count,
        )
        return False

    # Log the seq
    logger.info(
        "{newseq}\t{count}\t<green>ADDED</>",
        newseq=report_seq,
        count=seqstatus.count,
    )

    return True


def repair(
    config_path: pathlib.Path,
    msa_path: pathlib.Path,
    bedfile_path: pathlib.Path,
    output_dir: pathlib.Path,
    cores: int,
    force: bool,
):
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()  # Keep absolute path

    # Read in the config file
    with open(config_path) as f:
        base_cfg = json.load(f)

    # Overwrite the config file
    base_cfg["n_cores"] = cores  # Overwrite the number of cores
    base_cfg["output_dir"] = str(OUTPUT_DIR)
    base_cfg["primer_max_walk"] = 50
    base_cfg["primer_hairpin_th_max"] = 49.5

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Set up the logger
    logger = setup_loger(OUTPUT_DIR)

    # Read in the MSA file
    msa = MSA(
        name=msa_path.stem,
        path=msa_path,
        msa_index=0,
        mapping=base_cfg["mapping"],
        logger=logger,
    )
    if "/" in msa._chrom_name:
        new_chromname = msa._chrom_name.split("/")[0]
        logger.warning(
            "<red>WARNING</>: Having a '/' in the chromname {msachromname} will cause issues with figure generation bedfile output. Parsing chromname <yellow>{msachromname}</> -> <green>{new_chromname}</>",
            msachromname=msa._chrom_name,
            new_chromname=new_chromname,
        )
        msa._chrom_name = new_chromname

    logger.info(
        "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
        msa_path=msa.name,
        msa_rows=msa.array.shape[0],
        msa_cols=msa.array.shape[1],
    )

    # Update the base_cfg with the new msa
    # Create MSA checksum
    with open(msa_path, "rb") as f:
        msa_checksum = hashlib.file_digest(f, "md5").hexdigest()
    current_msa_index = max([int(x) for x in base_cfg["msa_data"].keys()])
    base_cfg["msa_data"][str(current_msa_index + 1)] = {
        "msa_name": msa.name,
        "msa_path": str("work/" + msa_path.name),
        "msa_chromname": msa._chrom_name,
        "msa_uuid": msa._uuid,
        "msa_checksum": msa_checksum,
    }
    # Copy the MSA file to the work dir
    local_msa_path = OUTPUT_DIR / "work" / msa_path.name
    shutil.copy(msa_path, local_msa_path)

    # Update the config with new algoversion
    base_cfg["algorithmversion"] = f"primalscheme3:{__version__}"

    # Read in the bedfile
    all_primerpairs, _header = read_in_bedprimerpairs(bedfile_path)

    # Get the primerpairs for this new MSA
    primerpairs_in_msa = [
        BedPrimerPairRepair(pp)
        for pp in all_primerpairs
        if pp.chrom_name == msa._chrom_name
    ]
    if len(primerpairs_in_msa) == 0:
        logger.error(
            "No primerpairs found for {msa._chrom_name} in {bedfile_path}",
            msa=msa,
            bedfile_path=bedfile_path,
        )
        sys.exit(1)

    # Get all the seqs in each pool
    seqs_in_pools = [[] for _ in range(base_cfg["npools"])]
    for pp in primerpairs_in_msa:
        seqs_in_pools[pp.pool].extend([*pp.fprimer.seqs, *pp.rprimer.seqs])

    # Find the indexes in the MSA that the primerbed refere to
    assert msa._mapping_array is not None
    mapping_list = list(msa._mapping_array)

    # For primerpair in the bedfile, check if new seqs need to be added by digestion the MSA
    for pp in primerpairs_in_msa:
        msa_fkmer_end = mapping_list.index(pp.fprimer.end)

        logger.info(
            "Checking <blue>{pp.amplicon_prefix}_{pp.amplicon_number}_LEFT</>",
            pp=pp,
        )
        _end_col, fseq_counts = f_digest_to_count(
            (msa.array, base_cfg, msa_fkmer_end, base_cfg["minbasefreq"])
        )
        # Check for early return conditions
        if detect_early_return(fseq_counts):
            logger.warning(
                "Early return for {pp.amplicon_prefix}_{pp.amplicon_number}_LEFT",
                pp=pp,
            )
            continue

        # Thermo check each sequence
        seqstatuss: list[SeqStatus] = []
        for seq, count in fseq_counts.items():
            if isinstance(seq, DIGESTION_ERROR):
                thermo_status = seq
                seq = None
            else:
                thermo_status = thermo_check_seq(seq, base_cfg)
            seqstatuss.append(SeqStatus(seq, count, thermo_status))

        # Decide if the new seqs should be added
        for seqstatus in seqstatuss:

            if not report_check(
                seqstatus=seqstatus,
                current_primer_seqs=pp.fprimer.seqs,
                seqs_in_pools=seqs_in_pools,
                pool=pp.pool,
                dimerscore=base_cfg["dimerscore"],
                logger=logger,
            ):
                continue

            # Add the new seq
            pp.new_fkmer_seqs.add(seqstatus.seq)
            seqs_in_pools[pp.pool].append(seqstatus.seq)

        # Handle the right primer
        logger.info(
            "Checking <blue>{pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT</>",
            pp=pp,
        )
        msa_rkmer_start = mapping_list.index(pp.rprimer.start)
        _start_col, rseq_counts = r_digest_to_count(
            (msa.array, base_cfg, msa_rkmer_start, base_cfg["minbasefreq"])
        )
        # Check for early return conditions
        if detect_early_return(rseq_counts):
            logger.warning(
                "Early return for {pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT",
                pp=pp,
            )
            continue

        # Valid seqs
        valid_rseqs = {
            reverse_complement(seq): count
            for seq, count in rseq_counts.items()
            if type(seq) == str
        }

        # Thermo check each sequence
        rseqstatuss: list[SeqStatus] = []
        for seq, count in valid_rseqs.items():
            if isinstance(seq, DIGESTION_ERROR):
                thermo_status = seq
                seq = None
            else:
                thermo_status = thermo_check_seq(seq, base_cfg)
            rseqstatuss.append(SeqStatus(seq, count, thermo_status))

        # Decide if the new seqs should be added
        for rseqstatus in rseqstatuss:
            if not report_check(
                seqstatus=rseqstatus,
                current_primer_seqs=pp.rprimer.seqs,
                seqs_in_pools=seqs_in_pools,
                pool=pp.pool,
                dimerscore=base_cfg["dimerscore"],
                logger=logger,
            ):
                continue

            # Add the new seq
            pp.new_rkmer_seqs.add(rseqstatus.seq)
            seqs_in_pools[pp.pool].append(rseqstatus.seq)

    # Write out the new bedfile
    with open(OUTPUT_DIR / "primer.bed", "w") as f:
        for pp in primerpairs_in_msa:
            f.write(pp.to_bed(new_primer_prefix=msa._uuid) + "\n")

    # Amplicon and primertrimmed files should not have changed. Can be copied from the input dir
    # Not sure how to handle the amplicon names, as the primerstem has changed?
    ## Keep orginal names for now

    # Write the config dict to file
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(base_cfg, sort_keys=True))
