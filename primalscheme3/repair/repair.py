import json
import pathlib
import sys

# Core imports
from primalscheme3.core.msa import MSA
from primalscheme3.core.bedfiles import read_in_bedprimerpairs, BedPrimerPair
from primalscheme3.core.logger import setup_loger
from primalscheme3.core.digestion import (
    f_digest_to_count,
    r_digest_to_count,
    DIGESTION_ERROR,
)
from primalscheme3.core.thermo import forms_hairpin, thermo_check_kmers, THERMORESULT
from primalscheme3.core.classes import FKmer, RKmer
from primalscheme3.core.seq_functions import reverse_complement

from primaldimer_py import do_pools_interact_py  # type: ignore


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


def repair(args):
    OUTPUT_DIR = pathlib.Path(args.output).absolute()  # Keep absolute path

    # Read in the config file
    with open(args.config) as f:
        base_cfg = json.load(f)

    # Overwrite the config file
    base_cfg["n_cores"] = args.cores  # Overwrite the number of cores
    base_cfg["output_dir"] = str(OUTPUT_DIR)
    base_cfg["primer_max_walk"] = 50
    base_cfg["primer_hairpin_th_max"] = 49.5

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Set up the logger
    logger = setup_loger(OUTPUT_DIR)

    # Read in the MSA file
    msa = MSA(
        name=args.msa.stem,
        path=args.msa,
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

    # Read in the bedfile
    all_primerpairs, _header = read_in_bedprimerpairs(args.bedfile)

    # Get the primerpairs for this new MSA
    primerpairs_in_msa = [
        BedPrimerPairRepair(pp)
        for pp in all_primerpairs
        if pp.chrom_name == msa._chrom_name
    ]
    if len(primerpairs_in_msa) == 0:
        logger.error(
            "No primerpairs found for {msa._chrom_name} in {args.bedfile}",
            msa=msa,
            args=args,
        )
        sys.exit(1)

    # Get all the seqs in each pool
    seqs_in_pools = [[] for _ in range(base_cfg["npools"])]
    for pp in primerpairs_in_msa:
        seqs_in_pools[pp.pool].extend([*pp.fprimer.seqs, *pp.rprimer.seqs])

    # Find the indexes in the MSA that the primerbed refere to
    assert msa._mapping_array is not None
    mapping_list = list(msa._mapping_array)

    # For pp in primerpairs_in_msa
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
        thermo_status = dict()
        for seq, count in fseq_counts.items():
            thermo_status[seq] = thermo_check_seq(seq, base_cfg)

        # Decide if the new seqs should be added
        for newseq, count in fseq_counts.items():
            # Check it passed thermo
            if thermo_status[newseq] != THERMORESULT.PASS:
                logger.warning(
                    "{newseq}\t<red>FAILED</>: {thermo_status[newseq]}",
                    thermo_status=thermo_status,
                    newseq=newseq.rjust(40),
                )
                continue

            # Check it is a new seq
            if newseq in pp.fprimer.seqs:
                logger.warning(
                    "{newseq}\t{count}\t<green>NOTADDED</>: Already Included",
                    newseq=newseq.rjust(40),
                    count=count,
                )
                continue

            # Check for minor allele
            if count < 5:
                logger.warning(
                    "{newseq}\t{count}\t<red>FAILED</>: Minor Allele",
                    newseq=newseq.rjust(40),
                    count=count,
                )
                continue

            # Check for dimer with pool
            if do_pools_interact_py(
                [newseq], seqs_in_pools[pp.pool], base_cfg["dimerscore"]
            ):
                logger.warning(
                    "{newseq}\t{count}\t<red>FAILED</>: Interaction with pool",
                    newseq=newseq.rjust(40),
                    count=count,
                )
                continue

            # Add the new seq
            pp.new_fkmer_seqs.add(newseq)
            seqs_in_pools[pp.pool].append(newseq)
            logger.info(
                "{newseq}\t{count}\t<green>ADDED</>",
                newseq=newseq.rjust(40),
                count=count,
            )
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
        thermo_status = dict()
        for seq, count in valid_rseqs.items():
            thermo_status[seq] = thermo_check_seq(seq, base_cfg)

        # Decide if the new seqs should be added
        for newseq, count in valid_rseqs.items():
            # Check it passed thermo
            if thermo_status[newseq] != THERMORESULT.PASS:
                logger.warning(
                    "{newseq}\t<red>FAILED</>: {thermo_status[newseq]}",
                    thermo_status=thermo_status,
                    newseq=newseq.rjust(40),
                )
                continue

            # Check it is a new seq
            if newseq in pp.rprimer.seqs:
                logger.warning(
                    "{newseq}\t{count}\t<green>NOTADDED</>: Already Included",
                    newseq=newseq.rjust(40),
                    count=count,
                )
                continue

            # Check for minor allele
            if count < 5:
                logger.warning(
                    "{newseq}\t{count}\t<red>FAILED</>: Minor Allele",
                    newseq=newseq.rjust(40),
                    count=count,
                )
                continue

            # Check for dimer with pool
            if do_pools_interact_py(
                [newseq], seqs_in_pools[pp.pool], base_cfg["dimerscore"]
            ):
                logger.warning(
                    "{newseq}\t{count}\t<red>FAILED</>: Interaction with pool",
                    newseq=newseq.rjust(40),
                    count=count,
                )
                continue

            # Add the new seq
            pp.new_rkmer_seqs.add(newseq)
            seqs_in_pools[pp.pool].append(newseq)
            logger.info(
                "{newseq}\t{count}\t<green>ADDED</>",
                newseq=newseq.rjust(40),
                count=count,
            )

    # Write out the new bedfile
    with open(OUTPUT_DIR / f"primer.bed", "w") as f:
        for pp in primerpairs_in_msa:
            f.write(pp.to_bed(new_primer_prefix=msa._uuid) + "\n")

    # Write the config dict to file
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        outfile.write(json.dumps(base_cfg, sort_keys=True))

    ## DO THIS LAST AS THIS CAN TAKE A LONG TIME
    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        OUTPUT_DIR / "work",
        last_pp_added=all_primerpairs,
    )
    generate_all_plots(plot_data, OUTPUT_DIR)
