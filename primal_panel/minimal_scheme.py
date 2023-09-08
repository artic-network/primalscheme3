# Module imports
from primal_panel.cli import cli
from primal_panel.minimal_scheme_classes import MSA, panel

# Submodule imports
from primal_digest.config import config_dict
from primal_digest.mismatches import MatchDB
from primal_digest.create_reports import generate_plot

# General import
import pathlib
import sys
from loguru import logger
import json


logger = logger.opt(colors=True)


def read_bedfile(path):
    ## If bedfile given, parse it:
    bed_lines = []
    try:
        with open(path, "r") as bedfile:
            for line in bedfile.readlines():
                line = line.strip()
                if line:
                    data = line.split()[:3]

                    # Validate bedfile
                    if int(data[1]) > int(data[2]):
                        raise ValueError(
                            f"Start position is greater than end position: {line}"
                        )
                    bed_lines.append(line.split()[:3])
    except:
        raise ValueError(f"Cannot parse bedfile: {path}")
    return bed_lines


def main():
    args = cli()
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict
    # Primer Digestion settings
    cfg["primer_gc_min"] = args.primer_gc_min
    cfg["primer_gc_max"] = args.primer_gc_max
    cfg["primer_tm_min"] = args.primer_tm_min
    cfg["primer_tm_max"] = args.primer_tm_max
    cfg["dimerscore"] = args.dimerscore
    cfg["n_cores"] = args.cores
    cfg["output_prefix"] = args.prefix
    cfg["output_dir"] = str(OUTPUT_DIR)
    cfg["msa_paths"] = [str(x) for x in ARG_MSA]
    cfg["amplicon_size_max"] = args.ampliconsizemax
    cfg["amplicon_size_min"] = args.ampliconsizemin
    cfg["min_overlap"] = args.minoverlap
    cfg["force"] = args.force
    cfg["npools"] = args.npools
    cfg["reducekmers"] = args.reducekmers
    cfg["minbasefreq"] = args.minbasefreq
    cfg["msa_index_to_ref_name"] = {
        index: msa_name for index, msa_name in enumerate(cfg["refname"])
    }
    # Add the mismatch params to the cfg
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 20
    cfg["mismatch_product_size"] = args.ampliconsizemax

    # Add plots to the cfg
    cfg["plot"] = True
    cfg["disable_progress_bar"] = False

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    ## Set up the loggers
    logger.remove()  # Remove default stderr logger
    # Add the deep log file
    logger.add(
        OUTPUT_DIR / "work/file.log",
        colorize=False,
        format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}",
        enqueue=True,
    )
    # Add the nice stdout progress
    logger.add(sys.stdout, colorize=True, format="{message}", level="INFO")

    # Create the mismatch db
    logger.info(
        "Creating the Mismatch Database",
    )
    mismatch_db = MatchDB(OUTPUT_DIR / "work/mismatch", ARG_MSA, cfg["primer_size_min"])
    logger.info(
        "<green>Created:</> {path}",
        path=f"{OUTPUT_DIR}/work/mismatch.db",
    )

    # Read in the bedfile if given
    if args.bedfile:
        msa_regions = {msa.name: set() for msa in ARG_MSA}
        bed_lines = read_bedfile(args.bedfile)

        # Map each bedfile to an MSA regions
        for bedline in bed_lines:
            # If the region can be mapped to an MSA
            if bedline[0] in msa_regions.keys():
                msa_regions.get(bedline[0]).add((int(bedline[1]), int(bedline[2])))
                logger.info(
                    "Mapped <green>{bedline_name}</>:<blue>{bedline_start}</>:<blue>{bedline_end}</> to <blue>{msa}</>",
                    bedline_name=bedline[0],
                    bedline_start=bedline[1],
                    bedline_end=bedline[2],
                    msa=bedline[0],
                )
            else:
                logger.info(
                    "Cannot map <red>{bedline_name}</>:<blue>{bedline_start}</>:<blue>{bedline_end}</>",
                    bedline_name=bedline[0],
                    bedline_start=bedline[1],
                    bedline_end=bedline[2],
                    msa=bedline[0],
                )

        # Turn the regions into a range of indexes
        msa_name_to_indexes = {}
        for msa_name, indexes in msa_regions.items():
            findexes = {
                fi
                for fi in (
                    range(x[0] - cfg["amplicon_size_max"], x[0] - 1) for x in indexes
                )
                for fi in fi
            }
            rindexes = {
                ri
                for ri in (
                    range(x[1] + 1, x[0] + cfg["amplicon_size_max"]) for x in indexes
                )
                for ri in ri
            }
            msa_name_to_indexes[msa_name] = (findexes, rindexes)

    ## Read in the MSAs
    msa_list: list[MSA] = []
    for msa_index, msa_path in enumerate(ARG_MSA):
        # Reads in the MSA
        msa = MSA(name=msa_path.name, path=msa_path, msa_index=msa_index)

        logger.info(
            "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
            msa_path=msa_path.name,
            msa_rows=msa.array.shape[0],
            msa_cols=msa.array.shape[1],
        )
        # Digest the MSA
        ## TODO Entry Point for targeted digestion
        if args.bedfile:
            msa.digest(cfg, indexes=msa_name_to_indexes.get(msa_path.name, False))
        else:
            msa.digest(cfg, indexes=False)
        logger.info(
            "<blue>{msa_path}</>: digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa_path.name,
            num_fkmers=len(msa.fkmers),
            num_rkmers=len(msa.rkmers),
        )
        # Generate all primerpairs
        msa.generate_primerpairs(cfg=cfg)
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa_path.name,
            num_pp=len(msa.primerpairs),
        )
        msa.primerpairs.sort(key=lambda x: -msa.get_pp_entropy(x) / len(x.all_seqs()))

        # Add the MSA to the list if it has primerpairs
        msa_list.append(msa)

    # Create a lookup dict for the msa index to name
    msa_index_to_name = {k: v for k, v in enumerate([msa.name for msa in ARG_MSA])}

    # Create the panel object
    panel = panel(msa_list, cfg, mismatch_db)
    counter = 0
    while panel.try_add_primerpair() and counter < len(ARG_MSA) * 3:
        logger.info(
            "Added <blue>amplicon</> for <blue>{msa_name}</>: {primer_start}\t{primer_end}",
            primer_start=panel._pool[-1].start,
            primer_end=panel._pool[-1].end,
            msa_name=msa_index_to_name.get(panel._pool[-1].msa_index),
        )
        counter += 1

    # Write primer bed file
    with open(OUTPUT_DIR / f"{cfg['output_prefix']}.primer.bed", "w") as outfile:
        primer_bed_str = []
        for pp in panel._pool:
            st = pp.__str__(
                msa_index_to_name.get(pp.msa_index, "NA"),
                msa_index_to_name.get(pp.msa_index, "NA"),
            )
            primer_bed_str.append(st.strip())
        outfile.write("\n".join(primer_bed_str))

    # Create all the plots
    if cfg["plot"]:
        for msa in panel.msas:
            generate_plot(msa, panel._pool, OUTPUT_DIR)

    # Write the config dict to file
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))


if __name__ == "__main__":
    main()
