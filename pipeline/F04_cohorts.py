##### BEGIN F04 cohorts, Table 1 and feasibility gate (BLIND) #####
# Builds the Q1 and Q2 analysed cohorts with protocol defaults. BLIND (STUDY-LIFECYCLE stage 7-8, P§12.1):
# events are shown in total only, never by exposure group, and no model is fitted. The outcome models (F05-F10)
# refuse to run until this output is saved in workbench_outputs/ and named in _config.RUN["blind_frozen_output"].
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import cohort as B
import outcomes as O

C.begin("F04", "cohorts, Table 1, feasibility (blind)")
NAMES = ["dx", "px", "drug", "visits", "labs", "persons", "death", "enroll", "survey", "ehr_span", "activity_months",
         "episodes", "courses", "prophylaxis", "egfr", "esrd", "ckd_lab", "rcc"]
D = {n: C.load_state(n) for n in NAMES}
D["atc_drugs"] = C.load_state("atc_drugs", internal=False)
G = C.load_state("drug_glossary", internal=False)
D["rcc_drug_starts"] = D["drug"].merge(G.loc[G["rcc_drug"], ["concept_id"]], left_on="drug_concept_id",
                                       right_on="concept_id")[["person_id", "start_date"]].rename(columns={"start_date": "date"})
SPECS = {"Q1": B.default_spec("Q1"), "Q2": B.default_spec("Q2")}
D_RCC_SENS = C.load_state("rcc_without_nephroureterectomy")
EXPO = {"Q1": "q1_group", "Q2": "q2_class"}


def _smd_bin(p1, p0):
    d = np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    return float((p1 - p0) / d) if d > 0 else 0.0


def table1(F, expo, levels, pairs=()):
    """Counts per level and group, plus standardized differences for each (group, reference) pair (covariates only)."""
    rows = []
    for v in K.TABLE1_VARS:
        if v not in F.columns:
            continue
        s = F[v]
        if pd.api.types.is_bool_dtype(s):
            s = s.map({True: "Yes", False: "No"})
        if pd.api.types.is_numeric_dtype(s) and s.nunique() > 8:
            r = {"var": v, "level": "p25 / p50 / p75"}
            for g in levels:
                r[f"n_{g}"] = np.nan
                r[f"pctl_{g}"] = C.pctl(F.loc[F[expo].eq(g), v])
            r["n_missing_all"] = int(F[v].isna().sum())
            for g1, g0 in pairs:
                x1, x0 = (pd.to_numeric(F.loc[F[expo].eq(g), v], errors="coerce").dropna() for g in (g1, g0))
                sd = np.sqrt((x1.var() + x0.var()) / 2) if len(x1) > 1 and len(x0) > 1 else np.nan
                r[f"smd_{g1}_vs_{g0}"] = float((x1.mean() - x0.mean()) / sd) if sd and sd > 0 else np.nan
            rows.append(r)
            continue
        sl = s.fillna("Missing").astype(str)
        for lev in sorted(sl.unique()):
            r = {"var": v, "level": lev}
            for g in levels:
                r[f"n_{g}"] = int((F[expo].eq(g) & sl.eq(lev)).sum())
            for g1, g0 in pairs:
                n1, n0 = int(F[expo].eq(g1).sum()), int(F[expo].eq(g0).sum())
                r[f"smd_{g1}_vs_{g0}"] = (_smd_bin(r[f"n_{g1}"] / n1, r[f"n_{g0}"] / n0) if n1 and n0 else np.nan)
            rows.append(r)
    return pd.DataFrame(rows)


def run_q(qn):
    def fn():
        F, flow = B.build(D, SPECS[qn])
        C.save_state(F, f"frame_{qn}")
        prev = C.CFG.get("footnote_extra", "")
        C.CFG["footnote_extra"] = (prev + " | " if prev else "") + F.attrs["spec"]
        C.show("flow", f"{qn} cohort flow (masked)", C.flow_mask(flow), counts=["n", "excluded"])
        levels = (["recurrent", "single", "intermediate"] if qn == "Q1" else ["A", "B", "C", "D"])
        T = table1(F, EXPO[qn], levels, K.SMD_PAIRS[qn])
        cnt = [c for c in T.columns if c.startswith("n_")]
        smd_deps = {f"smd_{g1}_vs_{g0}": [f"n_{g1}", f"n_{g0}"] for g1, g0 in K.SMD_PAIRS[qn]}
        C.show("table1", f"{qn} Table 1 by exposure group with standardized differences (masked; no outcome)", T,
               counts=cnt, group_cols=["var"], deps=smd_deps, max_rows=500,
               note="smd = standardized difference vs the reference group, before weighting (|smd| > 0.1 = imbalance); "
                    "blanked where a count is masked.")
        sizes = F[EXPO[qn]].value_counts().rename_axis("group").rename("n").reset_index()
        C.show("arms", f"{qn} analysed persons per exposure group (no outcome)", sizes, counts=["n"], group_cols=[])
        fu = pd.DataFrame([{"follow_up_years_p25_p50_p75": C.pctl(F["t_years"]),
                            "t0_year_p25_p50_p75": C.pctl(F["t0_year"]),
                            "years_t0_to_enrollment_p25_p50_p75": C.pctl(F["years_t0_to_enrollment"]),
                            "share_t0_before_enrollment": float((F["years_t0_to_enrollment"] > 0).mean())}])
        C.show("timing", f"{qn} follow-up and timing (all analysed persons; AL5, AL6)", fu, counts=[])
        # BLIND totals: events and person-years over the whole cohort, never by exposure
        # Sensitivity: rebuild with nephroureterectomy not confirming; report only how many events depend on it
        # (masked when 1-20, so the second total can never be derived).
        # Q1: totals over the compared groups only (recurrent + single), the same population as the flow's last step
        # and the gate. Totals that include the intermediate group would let its event count be derived by subtraction.
        F_s, _ = B.build(dict(D, rcc=D_RCC_SENS), SPECS[qn])
        if qn == "Q1":
            F = F.loc[F["q1_group"].isin(["recurrent", "single"])]
            F_s = F_s.loc[F_s["q1_group"].isin(["recurrent", "single"])]
        lost = int(F["event"].sum()) - int(F_s["event"].sum())
        tot = pd.DataFrame([{"n": len(F), "events": int(F["event"].sum()), "person_years": float(F["t_years"].sum()),
                             "deaths_competing": int(F["death_competing"].sum()),
                             "events_lost_without_nephroureterectomy": lost if lost >= 0 else np.nan}])
        C.show("events_total", f"{qn} RCC events in the analysed cohort, total only (blind)" + (" - recurrent + single" if qn == "Q1" else ""), tot,
               counts=["n", "events", "deaths_competing", "events_lost_without_nephroureterectomy"],
               note="events_lost_without_nephroureterectomy = confirmed events that need a nephroureterectomy to be confirmed.")
        C.CFG["footnote_extra"] = prev
    return fn


C.part("1 Q1 cohort (recurrent vs single UTI, landmark E1 + 365 d)", run_q("Q1"))
C.part("2 Q2 cohort (antibiotic class for the same episode, new users)", run_q("Q2"))


def p_gate():
    """D17 go / no-go and the events-per-parameter gate (STUDY-LIFECYCLE stage 8), from blind totals."""
    rows = []
    for qn, adj in (("Q1", K.ADJ_Q1), ("Q2", K.ADJ_Q2)):
        if not C.has_state(f"frame_{qn}"):
            continue
        F = C.load_state(f"frame_{qn}")
        if qn == "Q1":
            F = F.loc[F["q1_group"].isin(["recurrent", "single"])]
        ev = int(F["event"].sum())
        X, names = C.design(F, [c for c in adj if c in F.columns], cat_ref=None)
        k = len(names) + 1
        go = ("GO (>=300)" if ev >= K.DECISIONS["D17_go_events"] else
              "REDUCED (100-299): state minimum detectable HR" if ev >= K.DECISIONS["D17_reduced_events"] else
              "DESCRIPTIVE + Q2b only (<100)")
        smallest_arm = int(F[EXPO[qn]].value_counts().min())
        rows.append({"question": qn, "events": ev, "parameters": k, "epp": ev / k,
                     "D17_rule": go, "epp_ge_10": ev / k >= 10, "smallest_arm_above_20": smallest_arm > 20})
    C.show("gate", "Feasibility gate (D17 and events per parameter; blind totals)", pd.DataFrame(rows), counts=["events"],
           deps={"epp": ["events"]}, note="epp = events per parameter; blanked when events are masked.")
C.part("3 feasibility gate", p_gate)
O.stamp_f04(globals().get("CODE_HASH", ""))     # F05+ check that the frozen output came from this F04
C.end()
##### END F04 #####
