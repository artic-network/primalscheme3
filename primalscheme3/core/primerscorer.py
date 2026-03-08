import pathlib

import joblib
import pandas as pd

from primalscheme3.core.config import Config
from primalscheme3.core.seq_functions import complement
from primalscheme3.core.thermo import calc_binding


def find_mismatch(ref, query) -> tuple[int, str, str, str] | None:
    """
    This assume both seqs are in 5'-> 3' and are 3' anchored. With only a single SNP diff.
    Mutations are represented in 5'to 3' with the mutation in the second position.

    5' AAC
       |.|
    3' TAG
    AA/TA GA/CA

    mut_base: ref > query
        T>A

    """
    idx = None
    for i, (rb, qb) in enumerate(zip(ref[::-1], query[::-1], strict=False)):
        # if a single mutation is found get the idx
        if rb != qb:
            idx = i
            break

    # Check mutation is found
    if idx is None:
        return None

    dinu1 = (
        query[len(query) - idx - 2 : len(query) - idx]
        + "/"
        + complement(ref[len(ref) - idx - 2 : len(ref) - idx])
    )
    dinu2 = (
        complement(ref[len(ref) - idx])
        + complement(ref[len(ref) - idx - 1])
        + "/"
        + query[len(query) - idx]
        + query[len(query) - idx - 1]
    )

    mut_base = complement(dinu1[4]) + ">" + dinu1[1]
    return (idx, dinu1, dinu2, mut_base)


class PrimerScorer:
    def __init__(self, model_path: pathlib.Path):
        # Load data
        self.model = joblib.load(model_path)

    def _score(self, idx, dg, ddg, p3_dg, dinu1, dinu2, idx_prop, mut_base):
        """
        Predicts the Delta Ct cycle shift for a primer mismatch.

                Expected dictionary format:
                {
                    "dg": -4.1, "ddg": 2.0, "p3_dg": -3.5,
                    "idx": 15, "idx_prop": 0.75,
                    "mut_base": "A>G", "dinu1": "GC/CC", "dinu2": "AT/TG"
                }
        """
        # Convert dictionary to DataFrame and predict
        primer_features = {
            "dg": dg,
            "ddg": ddg,
            "p3_dg": p3_dg,
            "idx": idx,
            "idx_prop": idx_prop,
            "mut_base": mut_base,
            "dinu1": dinu1,
            "dinu2": dinu2,
        }
        df = pd.DataFrame([primer_features])
        prediction = self.model.predict(df)
        return round(float(prediction[0]), 3)

    def score_seq(self, ref: str, query: str, config: Config) -> None | float:
        """
        Generates  used in the model
        """
        #  Append '.' to allow terminal dinus
        aref = "." + ref + "."
        aquery = "." + query + "."

        # calculate dgs
        query_dg = calc_binding(query, query, config).dg
        ddg = calc_binding(query, ref, config).dg

        # find the mismatch
        r = find_mismatch(aref, aquery)
        if r is None:
            print(aref, aquery)
            return None
        idx, dinu1, dinu2, mut_base = r
        idx -= 1  # Account for the terminal "."

        idx_prop = idx / len(query)

        # Calc params
        p3_seq = query[len(query) - idx :]
        if p3_seq:
            p3_dg = calc_binding(p3_seq, p3_seq, config).dg
        else:
            p3_dg = 0

        return self._score(
            idx=idx,
            dg=query_dg,
            idx_prop=idx_prop,
            ddg=ddg,
            p3_dg=p3_dg,
            dinu1=dinu1,
            dinu2=dinu2,
            mut_base=mut_base,
        )


if __name__ == "__main__":
    print(
        find_mismatch("CTAAACGAACTTTAAAATCTGTGTGGCTG", "TCTAAACGAACTTTAAAAACTGTGTGGCTG")
    )
    ps = PrimerScorer(
        pathlib.Path(
            "/Users/kentcg/qpcr-single-mm/dinucleotide/q5_gam_surrogate_training_data.joblib"
        )
    )

    r = ps.score_seq(
        "CTAAACGAACTTTAAAATCTGTGTGGCTG", "TCTAAACGAACTTTAAAAACTGTGTGGCTG", Config()
    )
    print(r)

    exit()


# q  TCTAAACGAACTTTAAAAACTGTGTGGCTG
#     |||||||||||||||||||||||||||||
# ref CTAAACGAACTTTAAAATCTGTGTGGCTG


# q  5'  TCTAAACGAACTTTAAAAACTGTGTGGCTG  3'
#         |||||||||||||||||||||||||||||
# # refC  GATTTGCTTGAAATTTTAGACACACCGAC
