import json
import pathlib
import sys

# Core imports
from primalscheme3.core.msa import MSA
from primalscheme3.core.bedfiles import read_in_bedprimerpairs
from primalscheme3.core.logger import setup_loger
from primalscheme3.core.digestion import (
    f_digest_to_count,
    r_digest_to_count,
    process_seqs,
    DIGESTION_ERROR,
)
from primalscheme3.core.thermo import (
    forms_hairpin,
    thermo_check_kmers,
)
from primalscheme3.core.classes import FKmer, RKmer

from primaldimer_py import do_pools_interact_py  # type: ignore


class FKmerFixer:
    stop: int
    seqs: set[str]
    new_seqs: set[str]
    new_seqs_scores: dict[str, tuple[int, str | DIGESTION_ERROR]]

    def __init__(
        self,
        end: int,
        seqs: set[str],
        new_seqs: set[str],
        new_seqs_scores: dict[str, tuple[int, str | DIGESTION_ERROR]],
    ) -> None:
        self.seqs = seqs
        self.end = end
        self.new_seqs = new_seqs
        self.new_seqs_scores = new_seqs_scores

    def add_seq(self, seq: str) -> None:
        # Moves a seq from new_seqs to seq
        self.new_seqs.remove(seq)
        self.new_seqs_scores.pop(seq)

        self.seqs.add(seq)


class RKmerFixer:
    start: int
    seqs: set[str]
    new_seqs: set[str]
    new_seqs_scores: dict[str, tuple[int, str | DIGESTION_ERROR]]

    def __init__(
        self,
        start: int,
        seqs: set[str],
        new_seqs_scores: dict[str, tuple[int, str | DIGESTION_ERROR]],
    ) -> None:
        self.seqs = seqs
        self.start = start
        self.new_seqs = set()
        self.new_seqs_scores = new_seqs_scores

    def add_seq(self, seq: str) -> None:
        # Moves a seq from new_seqs to seq
        self.new_seqs.remove(seq)
        self.new_seqs_scores.pop(seq)

        self.seqs.add(seq)


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
        pp for pp in all_primerpairs if pp.chrom_name == msa._chrom_name
    ]
    if len(primerpairs_in_msa) == 0:
        raise ValueError(
            f"No primerpairs found for {msa._chrom_name} in {args.bedfile}"
        )

    # Find the indexes in the MSA that the primerbed refere to
    fkmer_ends = {x.fprimer.end for x in primerpairs_in_msa}
    rkmer_starts = {x.rprimer.start for x in primerpairs_in_msa}
    assert msa._mapping_array is not None
    mapping_list = list(msa._mapping_array)
    mapped_fkmer_ends = {mapping_list.index(fkmer_end) for fkmer_end in fkmer_ends}
    mapped_rkmer_starts = {
        mapping_list.index(rkmer_start) for rkmer_start in rkmer_starts
    }
    # Go lower into digestion process so we can extract the Seqs/Errors and there counts
    # For Fkmers

    repair_fkmers: list[FKmerFixer] = []
    for mapped_fkmer_end in mapped_fkmer_ends:
        errored = False
        _end_col, seq_counts = f_digest_to_count(
            (msa.array, base_cfg, mapped_fkmer_end, base_cfg["minbasefreq"])
        )
        # Check for early return conditions
        for error, count in seq_counts.items():
            if count == -1 and type(error) == DIGESTION_ERROR:
                errored = True
                print(f"errored: {error} {count}")
                break
        if errored:
            continue

        # Thermo check each sequence
        thermo_status = dict()
        for seq, count in seq_counts.items():
            # Thermo check
            if not thermo_check_kmers([seq], base_cfg):  # type: ignore
                thermo_status[seq] = DIGESTION_ERROR.THERMO_FAIL
                continue

            # Check for hairpins
            elif forms_hairpin([seq], base_cfg):  # type: ignore
                thermo_status[seq] = DIGESTION_ERROR.HAIRPIN_FAIL
                continue

            # Check for dimer
            if do_pools_interact_py([seq], [seq], base_cfg["dimerscore"]):
                thermo_status[seq] = DIGESTION_ERROR.DIMER_FAIL
                continue

            thermo_status[seq] = "PASS"

        print(_end_col)
        fkmer = [
            x.fprimer
            for x in primerpairs_in_msa
            if x.fprimer.end == msa._mapping_array[mapped_fkmer_end]
        ][0]
        fkmer_fixer = FKmerFixer(
            end=fkmer.end,
            seqs=fkmer.seqs,
            new_seqs={seq for seq in seq_counts.keys() if seq not in fkmer.seqs},  # type: ignore
            new_seqs_scores={seq: (count, thermo_status[seq]) for seq, count in seq_counts.items() if seq not in fkmer.seqs},  # type: ignore
        )

        print(len(fkmer_fixer.new_seqs_scores))
        print(len(fkmer_fixer.seqs))

        fkmer_fixer.add_seq([x for x in fkmer_fixer.new_seqs_scores.keys()][0])

        print(len(fkmer_fixer.new_seqs_scores))
        print(len(fkmer_fixer.seqs))

    exit()
    # Digest the MSA
    msa.digest(base_cfg)

    # Go lower to ensure digestion products are always returned even if error

    # Create a mapping array to map the msa to the digestion products

    # Get Fkmers that match the input scheme
    wanted_fkmer_ends = {x.fprimer.end for x in primerpairs_in_msa}
    new_fkmers = [fkmer for fkmer in msa.fkmers if fkmer.end in wanted_fkmer_ends]
    # Get Rkmers that match the input scheme
    wanted_rkmer_starts = {x.rprimer.start for x in primerpairs_in_msa}
    new_rkmers = [rkmer for rkmer in msa.rkmers if rkmer.start in wanted_rkmer_starts]

    print(len(new_fkmers), len(new_rkmers))
    print(new_fkmers[0].seqs)
    print(new_fkmers[0].end)
    exit()

    for x in primerpairs_in_msa:
        total_col_seqs = set()
        for row_index in range(0, msa.array.shape[0]):
            # Get the start and end col
            start_array = msa.array[
                row_index, x.fprimer.end - base_cfg["primer_size_min"] : x.fprimer.end
            ]
            start_seq = "".join(start_array).replace("-", "")

            results = wrap_walk(
                walk_left,
                msa.array,
                x.fprimer.end,
                x.fprimer.end - base_cfg["primer_size_min"],
                row_index=row_index,
                seq_str=start_seq,
                cfg=base_cfg,
            )
            total_col_seqs.update(results)
        # Thermo check each sequence

        seqs_results = {}
        for seq in total_col_seqs:
            if type(seq) != str:
                continue

            seqs_results[seq] = "PASS"

            # Check for hairpins
            if forms_hairpin([seq], base_cfg):
                seqs_results[seq] = "FAIL_HAIRPIN"
                continue
            # Check for THERMO
            if not thermo_check_kmers([seq], base_cfg):
                seqs_results[seq] = "FAIL_THERMO"
                continue
            # Check for dimer
            if do_pools_interact_py([seq], [seq], base_cfg["dimerscore"]):
                seqs_results[seq] = "FAIL_HOMODIMER"
                continue

        print(seqs_results)

    logger.info(
        "<blue>{msa_path}</>: digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
        msa_path=msa.name,
        num_fkmers=len(msa.fkmers),
        num_rkmers=len(msa.rkmers),
    )

    # Create a dodgy pool object to run the interaction checker on
    seqs_in_pools = [[] for _ in range(base_cfg["npools"])]
    for pp in all_primerpairs:
        seqs_in_pools[pp.pool].extend([*pp.fprimer.seqs, *pp.rprimer.seqs])

    # Keep track of what has been added
    addedseqs = set()
    # Repair the fprimer class in the primerpairs
    for pp in primerpairs_in_msa:
        # Get the fprimer from the msa for the same index
        newfprimer = [fkmer for fkmer in msa.fkmers if fkmer.end == pp.fprimer.end]
        if len(newfprimer) != 1:
            print(f"no new fprimer for {pp.fprimer.end}")
            continue

        # Find what new seqs need to be added
        new_fprimer_seqs = newfprimer[0].seqs.difference(pp.fprimer.seqs)
        # If new seqs run the interaction checker
        if new_fprimer_seqs:
            # guard for interaction
            if do_pools_interact_py(
                [*new_fprimer_seqs], seqs_in_pools[pp.pool], base_cfg["dimerscore"]
            ):
                logger.warning(
                    "<red>WARNING</>: Interaction detected for new fprimer {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                    pp_chromname=pp.chrom_name,
                    pp_ampliconprefix=pp.amplicon_prefix,
                    pp_amplicon_number=pp.amplicon_number,
                )
                raise ValueError()

            # Update the fprimer
            pp.fprimer.seqs.update(new_fprimer_seqs)
            addedseqs.update(new_fprimer_seqs)
            # Update the pool seqs
            seqs_in_pools[pp.pool].extend(new_fprimer_seqs)
            logger.info(
                "Added new seqs {seqs} to <green>fprimer</> {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                pp_chromname=pp.chrom_name,
                pp_ampliconprefix=pp.amplicon_prefix,
                pp_amplicon_number=pp.amplicon_number,
                seqs=",".join(new_fprimer_seqs),
            )

        # Find what rprimer seqs need to be added
        new_rprimer = [rkmer for rkmer in msa.rkmers if rkmer.start == pp.rprimer.start]
        if len(new_rprimer) != 1:
            raise ValueError(
                f"Could not find a new rprimer for {pp.chrom_name}:{pp.amplicon_prefix}:{pp.amplicon_number}:RIGHT"
            )
        new_rprimer_seqs = new_rprimer[0].seqs.difference(pp.rprimer.seqs)
        # If new seqs run the interaction checker
        if new_rprimer_seqs:
            # guard for interaction
            if do_pools_interact_py(
                [*new_rprimer_seqs], seqs_in_pools[pp.pool], base_cfg["dimerscore"]
            ):
                logger.warning(
                    "<red>WARNING</>: Interaction detected for new rprimer {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                    pp_chromname=pp.chrom_name,
                    pp_ampliconprefix=pp.amplicon_prefix,
                    pp_amplicon_number=pp.amplicon_number,
                )
                raise ValueError()

            # Update the rprimer
            pp.rprimer.seqs.update(new_rprimer_seqs)
            addedseqs.update(new_rprimer_seqs)
            # Update the pool seqs
            seqs_in_pools[pp.pool].extend(new_rprimer_seqs)
            logger.info(
                "Added new seqs {seqs} to <green>rprimer</> {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                pp_chromname=pp.chrom_name,
                pp_ampliconprefix=pp.amplicon_prefix,
                pp_amplicon_number=pp.amplicon_number,
                seqs=",".join(new_rprimer_seqs),
            )

    logger.info("Added {n} new seqs", n=len(addedseqs))
