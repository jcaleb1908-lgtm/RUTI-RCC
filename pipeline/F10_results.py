##### BEGIN F10 compiled results, robustness rule and glossary #####
# Collects every estimate F05-F09 saved (aggregate only), applies the P§12.1 robustness rule to the primary
# contrasts, and prints the variable glossary for the paper. Runs only after the freeze.
import numpy as np
import pandas as pd
import aou_core as C
import outcomes as O

C.begin("F10", "compiled results, robustness, glossary")
O.freeze_check()
RES = C.load_state("results", internal=False, required=False)
ROBUST_FROM = ("A9", "A10", "A11")      # P§12.1: lag / latency, outcome level, confirmation


def p_table():
    """All hazard ratios with E-values (point and CI limit closer to 1; rare outcome, HR ~ RR)."""
    if RES is None or not len(RES):
        C.show("all", "Compiled results", pd.DataFrame([{"status": "no results saved: run F05-F09 first"}]), counts=[])
        return
    cols = ["analysis_id", "analysis", "question", "contrast", "HR", "lo", "hi", "p", "E_value", "E_value_ci", "missing", "m", "status"]
    T = RES[[c for c in cols if c in RES.columns]].sort_values(["analysis_id", "question", "analysis"], kind="stable")
    C.show("all", "Compiled hazard ratios, every analysis", T, counts=[], max_rows=400)
C.part("1 compiled table", p_table)


def p_robust():
    """P§12.1: a primary finding is robust if the direction is unchanged and the HR changes by < 20% in the
    lag/latency, outcome-level and confirmation sensitivity analyses."""
    if RES is None or not len(RES):
        return
    prim = RES.loc[RES["analysis_id"].isin(["A5", "A5b"]) & RES["status"].eq("shown")]
    rows = []
    for _, p in prim.iterrows():
        s = RES.loc[RES["analysis_id"].isin(ROBUST_FROM) & RES["question"].eq(p["question"]) & RES["contrast"].eq(p["contrast"])]
        shown = s.loc[s["status"].eq("shown")]
        same_dir = np.sign(np.log(shown["HR"])) == np.sign(np.log(p["HR"]))
        change = (shown["HR"] / p["HR"] - 1).abs()
        fails = shown.loc[~same_dir | (change >= 0.2), "analysis"].tolist()
        rows.append({"question": p["question"], "contrast": p["contrast"], "HR_primary": p["HR"],
                     "sensitivity_analyses_shown": len(shown), "not_shown": int((~s["status"].eq("shown")).sum()),
                     "max_relative_change": float(change.max()) if len(change) else np.nan,
                     "robust": "yes" if len(shown) and not fails else ("no" if fails else "not assessable"),
                     "failing": "; ".join(fails)})
    C.show("robust", "Robustness rule (P§12.1)", pd.DataFrame(rows), counts=[],
           note="direction unchanged and |HR / HR_primary - 1| < 0.2 in A9, A10 and A11")
C.part("2 robustness rule", p_robust)


GLOSSARY = [
    ("q1_group", "Q1 exposure in [E1, E1 + 365 d): recurrent (>= 2 episodes within 180 d or >= 3 in 365 d), single (E1 only), intermediate (2 episodes > 180 d apart)", "P§2.1, 5.4"),
    ("q2_class", "Q2 index antibiotic class: A nitrofurantoin/fosfomycin (reference), B TMP-SMX/trimethoprim, C fluoroquinolone, D beta-lactam", "P§2.2, 6.1"),
    ("event", "Confirmed RCC (C64 / 189.0 plus a second code 30-365 d later, a confirming procedure, or a kidney-cancer drug)", "P§8"),
    ("t_years", "Years from outcome-counting start (T0 + 2 y) to the earliest of RCC, death, last EHR activity, ESRD, 2025-01-01", "P§7"),
    ("IPTW", "Stabilized inverse-probability-of-treatment weights from a (multinomial) logistic model on the adjustment set, truncated at the 1st/99th percentile", "P§2, 9.2; O15"),
    ("site", "EHR site of the last EHR visit before T0; Cox stratum", "P§9.2"),
    ("MI", "Chained equations with predictive mean matching for BMI, eGFR (Q2), education and income; outcome indicator and Nelson-Aalen hazard as predictors; m = 20 primary", "P§12.2 #8, D16"),
    ("CIF", "Aalen-Johansen cumulative incidence of RCC with death competing, IPTW-weighted", "P§12.2 #4"),
    ("E-value", "Minimum confounder association (risk-ratio scale) with exposure and outcome that could explain the HR away; HR treated as RR (rare outcome)", "P§10"),
    ("RERI", "Relative excess risk due to interaction from a joint-category IPTW Cox model, reference non-Hispanic White unexposed; delta-method CI", "P§11.3"),
]


def p_answerable():
    """What each analysis could answer: an estimate was shown, or it stays descriptive because events were too few
    (<= 20 in an arm or cell, O19), or the model was not estimable. D17: with < 100 events every estimate is exploratory."""
    if RES is None or not len(RES):
        return
    R_ = RES.assign(outcome=np.where(RES["status"].eq("shown"), "estimate shown (exploratory)",
                    np.where(RES["status"].astype(str).str.startswith("not estimable"), "not estimable",
                             "descriptive only: too few events")))
    T = R_.groupby(["analysis_id", "analysis", "question"], sort=False)["outcome"].agg(
        lambda s: s.value_counts().to_dict()).reset_index()
    T["estimates_shown"] = R_.groupby(["analysis_id", "analysis", "question"], sort=False)["status"].apply(
        lambda s: int(s.eq("shown").sum())).values
    T["estimates_total"] = R_.groupby(["analysis_id", "analysis", "question"], sort=False).size().values
    T["answer"] = np.where(T["estimates_shown"].eq(0), "DESCRIPTIVE ONLY (insufficient events)",
                           np.where(T["estimates_shown"].lt(T["estimates_total"]), "PARTLY (some contrasts descriptive)",
                                    "ESTIMATED (exploratory)"))
    C.show("answerable", "What each analysis could answer", T.drop(columns="outcome"), counts=[], max_rows=200,
           note="Exploratory throughout: D17 needs >= 100 events for a reduced and >= 300 for a confirmatory analysis.")
C.part("3 what could be answered", p_answerable)


def p_glossary():
    C.show("glossary", "Glossary of analysis terms", pd.DataFrame(GLOSSARY, columns=["term", "definition", "source"]), counts=[])
C.part("4 glossary", p_glossary)
C.end()
##### END F10 #####
