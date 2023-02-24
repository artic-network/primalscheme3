from jinja2 import Environment, PackageLoader, select_autoescape
from primal_digest.classes import PrimerPair
from typing import Union
from base64 import b64encode
from gzip import compress


# This function exports functions to convert text or binary data to a data URI readable by
# igv.js.
# ALL CREDIT: https://github.com/igvteam/igv-reports by https://github.com/jrobinso
def get_data_uri(data: Union[str, bytes]) -> str:

    """
    Return a data uri for the input, which can be either a string or byte array
    """

    if isinstance(data, str):
        data = compress(data.encode())
        mediatype = "data:application/gzip"
    else:
        if data[0] == 0x1F and data[1] == 0x8B:
            mediatype = "data:application/gzip"
        else:
            mediatype = "data:application:octet-stream"

    enc_str = b64encode(data)

    data_uri = mediatype + ";base64," + str(enc_str)[2:-1]
    return data_uri


def amplicon_primer_gff_rows(list_primer_pair: list[PrimerPair], cfg: dict) -> list[list[tuple[str, ...]]]:
    rows = [[] for _ in {x.pool for x in list_primer_pair}]
    ref_id = cfg.get("refname")

    # Amplicons
    for n,pp in enumerate(list_primer_pair):
        colors = ["#CEB992", "#73937E", "#585563", "#5B2E48", "#471323"]
        row = [
                ref_id,
                "primal-digest",
                "mRNA",
                pp.fprimer.end,
                pp.rprimer.start,  # GFF is closed
                ".",
                ".",
                ".",
                f"ID=AMPLICON_{pp.amplicon_number + 1};Color={colors[(pp.pool - 1) % (len(colors) - 1)]}",
            ]
        rows[pp.pool].append(tuple(map(str, row)))
            
        # FPrimers
        for fseq in pp.fprimer.seqs:
                row = [
                    ref_id,
                    "primal-digest",
                    "exon",
                    pp.fprimer.end - len(fseq),
                    pp.fprimer.end,
                    ".",
                    ".",
                    ".",
                    f"ID=test;Name=test;Parent=AMP{pp.amplicon_number + 1};Note=POOL_{pp.pool}",
                ]
                rows[pp.pool].append(tuple(map(str, row)))
        # RPrimers    
        for rseq in pp.rprimer.seqs:
                row = [
                    ref_id,
                    "primal-digest",
                    "exon",
                    pp.rprimer.start,
                    pp.rprimer.start + len(rseq),
                    ".",
                    ".",
                    ".",
                    f"ID=ID;Name=ID;Parent=AMP{n+1};Note=POOL_{pp.pool}",
                ]
                rows[pp.pool].append(tuple(map(str, row)))
        
        return rows

def amplicon_primer_gff(a_gff_rows: list[list[tuple[str, ...]]]) -> list[str]:
        """
        Return GFF3 formatted multiline strs for all amplicons and associated primers in
        the scheme.
        """
        gff_tracks = []
        for pool in a_gff_rows:
            pool_tsv_rows = ["\t".join(map(str, row)) for row in pool]
            gff_tracks.append("##gff-version 3\n" + "\n".join(pool_tsv_rows))
        return gff_tracks

def write_report(list_pp, OUTPUT_DIR, cfg, referance_fasta) -> None:
        """
        Write HTML report and associated JSON

        referance_fasta: The first Genome in the msa as .format("fasta")
        """
        # Data JSON
        data = {
            "primerTracks": [
                get_data_uri(track) for track in amplicon_primer_gff(amplicon_primer_gff_rows(list_pp, cfg))
            ],
        }

        # HTML Template
        jinja_env = Environment(
            loader=PackageLoader("primal_digest"), autoescape=select_autoescape()
        )
        template = jinja_env.get_template("report_template.html")
        rendered = template.render(data=data)
        
        with open(OUTPUT_DIR / f"{cfg.get('output_prefix')}_report.html", "w") as fh:
            fh.write(rendered)