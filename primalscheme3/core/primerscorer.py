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

    def _score(self, features_list: list[dict]):
        """
        Predicts the Delta Ct cycle shift for a list of primer mismatch features.
        """
        if not features_list:
            return []
        df = pd.DataFrame(features_list)
        predictions = self.model.predict(df)
        return [round(float(p), 3) for p in predictions]

    def _get_features(self, ref: str, query: str, config: Config) -> dict | None:
        """
        Extracts features for the model from a ref/query pair.
        """
        aref = "." + ref + "."
        aquery = "." + query + "."

        # calculate dgs
        query_dg = calc_binding(query, query, config).dg
        ddg = calc_binding(query, ref, config).dg

        # find the mismatch
        r = find_mismatch(aref, aquery)
        if r is None:
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

        return {
            "dg": query_dg,
            "ddg": ddg,
            "p3_dg": p3_dg,
            "idx": idx,
            "idx_prop": idx_prop,
            "mut_base": mut_base,
            "dinu1": dinu1,
            "dinu2": dinu2,
        }

    def score_seq_batch(
        self, pairs: list[tuple[str, str]], config: Config
    ) -> list[float | None]:
        """
        Predicts Delta Ct for a list of (ref, query) pairs.
        """
        features_list = []
        valid_indices = []

        for i, (ref, query) in enumerate(pairs):
            features = self._get_features(ref, query, config)
            if features:
                features_list.append(features)
                valid_indices.append(i)

        if not features_list:
            return [None] * len(pairs)

        # Batch predict
        predictions = self._score(features_list)

        # Map back to original list
        results = [None] * len(pairs)
        for idx, pred in zip(valid_indices, predictions, strict=True):
            results[idx] = pred

        return results

    def score_seq(self, ref: str, query: str, config: Config) -> None | float:
        """
        Generates prediction for a single sequence.
        """
        features = self._get_features(ref, query, config)
        if features is None:
            return None

        results = self._score([features])
        return results[0] if results else None


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
