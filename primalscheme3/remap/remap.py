# These functions are used to map a primerscheme to a new reference genome

import pathlib
import sys

from Bio import Seq, SeqIO, SeqRecord

from primalscheme3.core.bedfiles import (
    read_in_bedlines,
)
from primalscheme3.core.logger import setup_logger


def remap(
    bedfile_path: pathlib.Path,
    msa_path: pathlib.Path,
    id_to_remap_to: str,
    output_dir: pathlib.Path,
):
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()  # Keep absolute path

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Setup logger
    logger = setup_logger(OUTPUT_DIR)

    # Read in the primerpairs
    bedlines, _header = read_in_bedlines(bedfile_path)

    # Read in the MSA
    msa_dict = SeqIO.index(str(msa_path), "fasta")

    # Check primer's reference genome is in the msa
    chrom_names = set([primer.chrom_name for primer in bedlines])
    if len(chrom_names) > 1:
        # TODO specify the primer chrom name to remap
        logger.error("More than one reference genome in the primer.bed")
        sys.exit(1)
    remap_chrom = next(iter(chrom_names))
    if remap_chrom not in msa_dict:
        logger.error(
            "<red>{remap_chrom}</> ID from primer.bed is <red>not found</> in the MSA",
            remap_chrom=remap_chrom,
        )
        sys.exit(1)
    else:
        logger.info(
            "<green>{remap_chrom}</> ID from primer.bed is <green>found</> in the MSA",
            remap_chrom=remap_chrom,
        )

    # Check the new reference genome is in the msa
    if id_to_remap_to not in msa_dict:
        logger.error(
            "<red>{id_to_remap_to}</> remapping ID is <red>not found</> in the MSA",
            id_to_remap_to=id_to_remap_to,
        )
        sys.exit(1)
    else:
        logger.info(
            "<green>{id_to_remap_to}</> remapping ID is <green>found</> in the MSA",
            id_to_remap_to=id_to_remap_to,
        )

    # The the primer's reference genome to the MSA index
    primer_to_msa: dict[int, int] = {}
    ref_index = 0
    for msa_index, ref_base in enumerate(msa_dict[remap_chrom]):  # type: ignore
        if ref_base not in {"", "-"}:
            primer_to_msa[ref_index] = msa_index
            ref_index += 1

    # Create a dict that can map MSA indexes to the new reference genome
    msa_to_new_ref: dict[int, int] = {}
    new_index = 0
    for msa_index, ref_base in enumerate(msa_dict[id_to_remap_to]):  # type: ignore
        if ref_base not in {"", "-"}:
            msa_to_new_ref[msa_index] = new_index
            new_index += 1

    # Grab genome length
    msa_length = len(msa_dict.get(id_to_remap_to))  # type: ignore
    perfect_map = True

    for bedline in bedlines:
        primername = "_".join(bedline.primername.split("_")[:3])

        if bedline.direction == "+":
            # Dict will always have the key
            pp_fp_msa = primer_to_msa[bedline.end]
            pp_fp_newref = msa_to_new_ref.get(pp_fp_msa)  # type: ignore

            if pp_fp_newref is None:
                logger.warning(
                    f"<red>Gap preventing direct mapping</> of {primername}_LEFT {remap_chrom}:{bedline.end} -> {id_to_remap_to}"
                )
                # Walk to next valid index
                while pp_fp_newref is None and pp_fp_msa < msa_length:
                    pp_fp_msa += 1
                    pp_fp_newref = msa_to_new_ref.get(pp_fp_msa)

                # Check fixed
                if pp_fp_newref is None:
                    logger.critical(
                        f"Could not find or repair a valid index for {primername}_LEFT: {bedline.end}"
                    )
                    sys.exit(1)
                else:
                    logger.warning(
                        f"Fixed with non-direct mapping {remap_chrom}:{bedline.end} -> {id_to_remap_to}:{pp_fp_newref}"
                    )
                    perfect_map = False

            logger.debug(
                f"<green>Mapped {primername}_LEFT</>: {remap_chrom}:{bedline.end} -> {id_to_remap_to}:{pp_fp_newref}"
            )

            bedline._end = pp_fp_newref
            bedline._start = bedline.end - len(bedline.sequence)
            bedline.chrom_name = id_to_remap_to
        else:
            # Map the reverse primer
            pp_rp_msa = primer_to_msa[bedline.start]
            pp_rp_newref = msa_to_new_ref.get(pp_rp_msa)  # type: ignore

            if pp_rp_newref is None:
                logger.warning(
                    f"<red>Gap preventing direct mapping</> of {primername}_RIGHT {remap_chrom}:{bedline.start} -> {id_to_remap_to}"
                )
                # Walk left to next valid index
                while pp_rp_newref is None and pp_rp_msa > 0:
                    pp_rp_msa -= 1
                    pp_rp_newref = msa_to_new_ref.get(pp_rp_msa)

                if pp_rp_newref is None:
                    logger.critical(
                        f"Could not find or repair a valid index for {primername}_RIGHT: {bedline.start}"
                    )
                    sys.exit(1)
                else:
                    logger.warning(
                        f"Fixed with non-direct mapping {remap_chrom}:{bedline.start} -> {id_to_remap_to}:{pp_rp_newref}"
                    )
                    perfect_map = False

            logger.debug(
                f"<green>Mapped {primername}_RIGHT</>: {remap_chrom}:{bedline.start} -> {id_to_remap_to}:{pp_rp_newref}"
            )
            bedline._start = pp_rp_newref
            bedline._end = bedline.start + len(bedline.sequence)
            bedline.chrom_name = id_to_remap_to

    # Check if the mapping was perfect
    if perfect_map:
        logger.info("<green>Exact</> mapping")
    else:
        logger.warning("<red>Imperfect</> mapping. See .log for details.")

    # Write out the output files
    _header.append(f"# remapped {remap_chrom} -> {id_to_remap_to}")
    _header.append(f"# exact_mapping: {perfect_map}")

    # Write the new primer.bed file
    with open(OUTPUT_DIR / "primer.bed", "w") as f:
        for bedline in bedlines:
            f.write(str(bedline) + "\n")

    # Write the new reference genome out
    with open(OUTPUT_DIR / "reference.fasta", "w") as f:
        records = [
            SeqRecord.SeqRecord(
                Seq.Seq(
                    str(msa_dict[id_to_remap_to].seq).strip().replace("-", "").upper(),  # type: ignore
                ),
                id=id_to_remap_to,
                description="",
            )
        ]
        SeqIO.write(
            records,
            f,
            "fasta",
        )

    logger.info("Writing output files...")
    logger.info("<green>Run complete</>")
