##### BEGIN F07 missing data: multiple imputation vs complete case vs missing category #####
# Protocol 12.2 analysis 8 and D16 (ANALYSIS-MAP A8). Primary (MI m = 20) is F05.cox; this cell adds the
# sensitivity versions and the imputation diagnostics. Runs only after the freeze.
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import outcomes as O

C.begin("F07", "missing data (A8)")
O.freeze_check()
FR = {q: C.load_state(f"frame_{q}") for q in ("Q1", "Q2")}


def p_compare():
    """A8: the primary contrasts under complete case and missing-as-category."""
    rows = []
    for q, aid in (("Q1", "A8"), ("Q2", "A8")):
        for method in ("complete_case", "category"):
            rows += O.contrast(FR[q], q, aid, f"{q} primary, {method.replace('_', ' ')}", method=method)
    O.show_rows("compare", "Missing-data sensitivity: complete case and missing category (A8)", rows,
                note="Primary multiple-imputation estimates are in F05.cox")
C.part("1 complete case and category (A8)", p_compare)


def p_diag():
    """Imputation diagnostics: missing counts and observed vs imputed means (first 5 imputations pooled)."""
    rows = []
    for q in ("Q1", "Q2"):
        F = FR[q].loc[FR[q][O.EXPO[q]].isin(O.LEVELS[q])].reset_index(drop=True)
        sets = O.impute(F, q, O.EXPO[q], 5)
        for c in O.targets(q):
            miss = O.is_missing(F, c).values
            enc = (lambda s: pd.to_numeric(s, errors="coerce")) if c in ("bmi", "egfr") else \
                  (lambda s: s.map({lev: i for i, lev in enumerate(K.ORDINAL_LEVELS[c])}))
            obs = enc(F.loc[~miss, c]).astype(float)
            imp = pd.concat([enc(d.loc[miss, c]).astype(float) for d in sets]) if miss.any() else pd.Series(dtype=float)
            rows.append({"question": q, "variable": c, "n_missing": int(miss.sum()), "share_missing": float(miss.mean()),
                         "mean_observed": float(obs.mean()) if len(obs) else np.nan,
                         "mean_imputed": float(imp.mean()) if len(imp) else np.nan})
    C.show("diag", "Imputation diagnostics (education and income on their ordinal code 0-4)", pd.DataFrame(rows),
           counts=["n_missing"], deps={"share_missing": ["n_missing"], "mean_imputed": ["n_missing"]},
           note="PMM with 5 donors; predictors: exposure, adjustment set, event indicator, Nelson-Aalen hazard")
C.part("2 imputation diagnostics", p_diag)
C.end()
##### END F07 #####
