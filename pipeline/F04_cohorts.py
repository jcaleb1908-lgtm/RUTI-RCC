##### BEGIN F04 cohorts, Table 1 and feasibility gate (BLIND) #####
# Builds the Q1 and Q2 analysed cohorts with protocol defaults. BLIND (STUDY-LIFECYCLE stage 7-8, P§12.1):
# events are shown in total only, never by exposure group, and no model is fitted. Outcome models (F05+) are
# written only after the PI signs the protocol and this output is frozen.
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import cohort as B

C.begin("F04", "cohorts, Table 1, feasibility (blind)")
NAMES = ["dx", "px", "drug", "visits", "labs", "persons", "death", "enroll", "survey", "ehr_span", "activity_months",
         "episodes", "courses", "prophylaxis", "egfr", "esrd", "ckd_lab", "rcc"]
D = {n: C.load_state(n) for n in NAMES}
D["atc_drugs"] = C.load_state("atc_drugs", internal=False)
G = C.load_state("drug_glossary", internal=False)
D["rcc_drug_starts"] = D["drug"].merge(G.loc[G["rcc_drug"], ["concept_id"]], left_on="drug_concept_id",
                                       right_on="concept_id")[["person_id", "start_date"]].rename(columns={"start_date": "date"})
SPECS = {"Q1": B.default_spec("Q1"), "Q2": B.default_spec("Q2")}
EXPO = {"Q1": "q1_group", "Q2": "q2_class"}


def table1(F, expo, levels):
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
            rows.append(r)
            continue
        for lev in sorted(s.fillna("Missing").astype(str).unique()):
            r = {"var": v, "level": lev}
            for g in levels:
                r[f"n_{g}"] = int((F[expo].eq(g) & s.fillna("Missing").astype(str).eq(lev)).sum())
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
        T = table1(F, EXPO[qn], levels)
        cnt = [c for c in T.columns if c.startswith("n_")]
        C.show("table1", f"{qn} Table 1 by exposure group (masked; no outcome)", T, counts=cnt, group_cols=["var"],
               max_rows=400)
        sizes = F[EXPO[qn]].value_counts().rename_axis("group").rename("n").reset_index()
        C.show("arms", f"{qn} analysed persons per exposure group (no outcome)", sizes, counts=["n"], group_cols=[])
        fu = pd.DataFrame([{"follow_up_years_p25_p50_p75": C.pctl(F["t_years"]),
                            "t0_year_p25_p50_p75": C.pctl(F["t0_year"]),
                            "years_t0_to_enrollment_p25_p50_p75": C.pctl(F["years_t0_to_enrollment"]),
                            "share_t0_before_enrollment": float((F["years_t0_to_enrollment"] > 0).mean())}])
        C.show("timing", f"{qn} follow-up and timing (all analysed persons; AL5, AL6)", fu, counts=[])
        # BLIND totals: events and person-years over the whole cohort, never by exposure
        tot = pd.DataFrame([{"n": len(F), "events": int(F["event"].sum()), "person_years": float(F["t_years"].sum()),
                             "deaths_competing": int(F["death_competing"].sum())}])
        C.show("events_total", f"{qn} RCC events in the analysed cohort, total only (blind)", tot,
               counts=["n", "events", "deaths_competing"])
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
C.end()
##### END F04 #####
