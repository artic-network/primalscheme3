import pathlib

import joblib
import pandas as pd

from primalscheme3.core.config import Config
from primalscheme3.core.seq_functions import complement
from primalscheme3.core.thermo import calc_binding

# Allawi/SantaLucia internal mismatch NN ΔG₃₇ (cal/mol)
# Computed as dH - 310.15 * dS
# Keys: XY_WZ in antiparallel orientation
INTERNAL_MM_DG = {
    "AA_AT": 3293,
    "AA_CT": 1320,
    "AA_GT": 706,
    "AA_TA": 673,
    "AA_TC": 873,
    "AA_TG": 113,
    "AC_AG": 140,
    "AC_CG": 478,
    "AC_GG": -492,
    "AC_TA": 773,
    "AC_TC": 1364,
    "AC_TT": 638,
    "AG_AC": 402,
    "AG_CC": 786,
    "AG_GC": 94,
    "AG_TA": 13,
    "AG_TG": -154,
    "AG_TT": 721,
    "AT_AA": 673,
    "AT_CA": 773,
    "AT_GA": 13,
    "AT_TC": 723,
    "AT_TG": 73,
    "AT_TT": 650,
    "CA_AT": 920,
    "CA_CT": 1014,
    "CA_GA": 402,
    "CA_GC": 753,
    "CA_GG": 13,
    "CA_TT": 783,
    "CC_AG": 797,
    "CC_CG": 841,
    "CC_GA": 786,
    "CC_GC": 733,
    "CC_GT": 595,
    "CC_TG": 1014,
    "CG_AC": 753,
    "CG_CC": 733,
    "CG_GA": 94,
    "CG_GG": -156,
    "CG_GT": -472,
    "CG_TC": 392,
    "CT_AA": 873,
    "CT_CA": 1364,
    "CT_GC": 392,
    "CT_GG": -282,
    "CT_GT": -102,
    "CT_TA": 723,
    "GA_AT": 483,
    "GA_CA": 140,
    "GA_CC": 797,
    "GA_CG": -290,
    "GA_GT": 484,
    "GA_TT": 343,
    "GC_AG": -290,
    "GC_CA": 478,
    "GC_CC": 841,
    "GC_CT": 625,
    "GC_GG": -1101,
    "GC_TG": -585,
    "GG_AC": 13,
    "GG_CA": -492,
    "GG_CG": -1101,
    "GG_CT": 75,
    "GG_GC": -156,
    "GG_TC": -282,
    "GG_TT": 645,
    "GT_AA": 113,
    "GT_CC": 1014,
    "GT_CG": -585,
    "GT_CT": -596,
    "GT_GA": -154,
    "GT_TA": 73,
    "GT_TG": 1153,
    "TA_AA": 3293,
    "TA_AC": 920,
    "TA_AG": 483,
    "TA_CT": 983,
    "TA_GT": 427,
    "TA_TT": 665,
    "TC_AA": 1320,
    "TC_AC": 1014,
    "TC_AT": 983,
    "TC_CG": 625,
    "TC_GG": 75,
    "TC_TG": -596,
    "TG_AA": 706,
    "TG_AG": 484,
    "TG_AT": 427,
    "TG_CC": 595,
    "TG_GC": -472,
    "TG_GT": 523,
    "TG_TC": -102,
    "TT_AC": 783,
    "TT_AG": 343,
    "TT_AT": 665,
    "TT_CA": 638,
    "TT_GA": 721,
    "TT_GG": 645,
    "TT_TA": 650,
}

# Terminal mismatch NN ΔG₃₇ (cal/mol)
TERMINAL_MM_DG = {
    "AA_AT": -547,
    "AA_CT": -529,
    "AA_GT": -602,
    "AA_TA": -682,
    "AA_TC": -360,
    "AA_TG": -536,
    "AC_AG": -1022,
    "AC_CG": -998,
    "AC_GG": -1064,
    "AC_TA": -621,
    "AC_TC": -55,
    "AC_TT": -373,
    "AG_AC": -982,
    "AG_CC": -839,
    "AG_GC": -1193,
    "AG_TA": -670,
    "AG_TG": -449,
    "AG_TT": -502,
    "AT_AA": -682,
    "AT_CA": -621,
    "AT_GA": -670,
    "AT_TC": -346,
    "AT_TG": -586,
    "AT_TT": -384,
    "CA_AT": -470,
    "CA_CT": -297,
    "CA_GA": -982,
    "CA_GC": -770,
    "CA_GG": -923,
    "CA_TT": -328,
    "CC_AG": -719,
    "CC_CG": -613,
    "CC_GA": -839,
    "CC_GC": -518,
    "CC_GT": -720,
    "CC_TG": -581,
    "CG_AC": -770,
    "CG_CC": -518,
    "CG_GA": -1193,
    "CG_GG": -854,
    "CG_GT": -1010,
    "CG_TC": -613,
    "CT_AA": -360,
    "CT_CA": -55,
    "CT_GC": -613,
    "CT_GG": -800,
    "CT_GT": -859,
    "CT_TA": -346,
    "GA_AT": -542,
    "GA_CA": -1022,
    "GA_CC": -719,
    "GA_CG": -857,
    "GA_GT": -263,
    "GA_TT": -561,
    "GC_AG": -857,
    "GC_CA": -998,
    "GC_CC": -613,
    "GC_CT": -712,
    "GC_GG": 5256,
    "GC_TG": -906,
    "GG_AC": -923,
    "GG_CA": -1064,
    "GG_CG": 5256,
    "GG_CT": -769,
    "GG_GC": -854,
    "GG_TC": -800,
    "GG_TT": 0,
    "GT_AA": -536,
    "GT_CC": -581,
    "GT_CG": -906,
    "GT_CT": -824,
    "GT_GA": -449,
    "GT_TA": -586,
    "GT_TG": 0,
    "TA_AA": -547,
    "TA_AC": -470,
    "TA_AG": -542,
    "TA_AT": -796,
    "TA_CT": -547,
    "TA_GT": -643,
    "TA_TT": -441,
    "TC_AA": -529,
    "TC_AC": -297,
    "TC_AT": -547,
    "TC_CG": -712,
    "TC_GG": -769,
    "TC_TG": -824,
    "TG_AA": -602,
    "TG_AG": 863,
    "TG_AT": -643,
    "TG_CC": -720,
    "TG_GC": -1010,
    "TG_TC": -859,
    "TT_AC": -328,
    "TT_AG": -561,
    "TT_AT": -441,
    "TT_CA": -373,
    "TT_GA": -502,
    "TT_TA": -384,
}

GC_RICH_DINU2 = {"CC/GC", "GC/CT", "GC/CC", "CT/GC"}


def find_mismatch(ref, query) -> tuple[int, str, str, str] | None:
    mismatches = []
    for i, (rb, qb) in enumerate(zip(ref[::-1], query[::-1], strict=False)):
        # ref/query are padded with "." at both ends. A position where one side
        # is padding and the other is a base is a 5'-length difference between two
        # 3'-anchored primers, NOT a substitution, so it must not be counted.
        # (Padding-vs-padding at the shared 3' terminus is a match anyway.)
        if rb == "." or qb == ".":
            continue
        if rb != qb:
            mismatches.append(i)

    # Check exactly one mutation is found
    if len(mismatches) != 1:
        return None

    idx = mismatches[0]

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


def lookup_mm_dg(dinu_key: str, is_terminal: bool) -> float:
    """Look up mismatch NN ΔG₃₇. Returns 0 if not found."""
    key = dinu_key.replace("/", "_").replace(".", "")
    table = TERMINAL_MM_DG if is_terminal else INTERNAL_MM_DG
    return table.get(key, 0)


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

        r = find_mismatch(aref, aquery)
        if r is None:
            return None

        idx, dinu1, dinu2, mut_base = r
        idx -= 1  # Account for terminal "."
        is_terminal = idx == 0
        idx_prop = idx / len(query)

        # p3_dg: stability of 3' segment
        p3_seq = query[len(query) - idx :]
        p3_dg = calc_binding(p3_seq, p3_seq, config).dg if p3_seq else 0

        # prdg: proofread correction energy
        # Correct the mutant base back to the reference base
        ref_base = complement(dinu1[4])  # original primer base
        proof_seq_l = list(query)
        proof_seq_l[len(query) - idx - 1] = ref_base
        proof_seq = "".join(proof_seq_l)
        prdg = calc_binding(query, proof_seq, config).dg

        # mm_dg2: mismatch NN ΔG at the 3' flanking context
        if is_terminal:
            mm_dg2 = 0
        else:
            mm_dg2 = lookup_mm_dg(dinu2, False)

        # dinu2_gc_rich
        dinu2_gc_rich = 1 if dinu2 in GC_RICH_DINU2 else 0

        return {
            "p3_dg": p3_dg,
            "prdg": prdg,
            "idx_prop": idx_prop,
            "mm_dg2": mm_dg2,
            "mut_base": mut_base,
            "dinu2_gc_rich": dinu2_gc_rich,
        }

    def score_seq_batch(
        self, pairs: list[tuple[str, str]], config: Config
    ) -> list[float | None]:
        features_list = []
        valid_indices = []
        results: list[float | None] = [None] * len(pairs)

        for i, (ref, query) in enumerate(pairs):
            if ref == query:
                results[i] = 0.0
                continue
            features = self._get_features(ref, query, config)
            if features:
                features_list.append(features)
                valid_indices.append(i)

        if not features_list:
            return results

        predictions = self._score(features_list)

        for idx, pred in zip(valid_indices, predictions, strict=True):
            results[idx] = pred

        return results

    def score_seq(self, ref: str, query: str, config: Config) -> None | float:
        if ref == query:
            return 0.0
        features = self._get_features(ref, query, config)
        if features is None:
            return None
        results = self._score([features])
        return results[0] if results else None
