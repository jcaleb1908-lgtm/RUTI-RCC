##### BEGIN F05 main estimates #####
# Protocol 12.2 analyses 3-7 (ANALYSIS-MAP A3-A7c). Runs only after the blind F00-F04 output is frozen.
# Primary: IPTW cause-specific Cox with multiple imputation (m = 20), EHR site as stratum, robust variance;
# Aalen-Johansen cumulative incidence with death competing. Masking: event counts by arm only in the rates table.
import pandas as pd
import aou_core as C
import _config as K
import outcomes as O

C.begin("F05", "main estimates (rates, cumulative incidence, IPTW Cox, PH, Q1 types, Q1b, head-to-head)")
O.freeze_check()
FR = {q: C.load_state(f"frame_{q}") for q in ("Q1", "Q2")}


def p_rates():
    """A3: events, person-years and rates per 1,000 person-years by exposure group."""
    for q, levels in (("Q1", ["recurrent", "single", "intermediate"]), ("Q2", ["A", "B", "C", "D"])):
        C.show(f"rates_{q}", f"{q} RCC events, person-years and rate per 1,000 person-years by group", O.rates(FR[q], q, levels=levels),
               counts=["n", "events"], group_cols=[], deps={"rate_per_1000py": ["events"]},
               note="person-years from outcome-counting start (T0 + 2 y). Complementary cells hidden across groups.")
C.part("1 rates (A3)", p_rates)


def p_cif():
    """A4: IPTW Aalen-Johansen cumulative incidence at 5 and 10 years after outcome counting starts, and the difference."""
    for q in ("Q1", "Q2"):
        C.show(f"cif_{q}", f"{q} cumulative incidence of RCC (%), death competing, IPTW, m = {K.ANALYSIS['mice_m']}",
               O.cif(FR[q], q, method="mice"), counts=[],
               note="CIF and RD in %; 95% CI from Rubin's rules with bootstrap within-imputation variance (weights fixed); "
                    "shown only when each arm has >= 21 at risk and >= 21 events by that time.")
C.part("2 cumulative incidence (A4)", p_cif)


def p_cox():
    """A5 / A5b: primary IPTW cause-specific Cox, multiple imputation m = 20 (D16)."""
    rows = O.contrast(FR["Q1"], "Q1", "A5", "Q1 primary: recurrent vs single UTI", method="mice")
    rows += O.contrast(FR["Q2"], "Q2", "A5b", "Q2 primary: antibiotic class vs nitrofurantoin/fosfomycin", method="mice")
    O.show_rows("cox", "Primary hazard ratios (IPTW cause-specific Cox, MI m = 20)", rows)
C.part("3 primary Cox (A5, A5b)", p_cox)


def p_ph():
    """A6: proportional hazards: score test on log time per exposure term; period HRs 2-5 and > 5 years after T0."""
    T = pd.concat([O.ph_check(FR[q], q) for q in ("Q1", "Q2")], ignore_index=True)
    C.show("ph", "Proportional-hazards check (missing-category data, IPTW)", T, counts=[],
           note="p_schoenfeld_logt: approximate score test of beta(t) = beta + theta*log(t); periods are years after T0 "
                "(outcome counting starts at 2 y). If p < 0.05 the period HRs are the ones to report (P§12.2 #6).")
C.part("4 PH check (A6)", p_ph)


def p_secondary():
    """A7a Q1 type-specific (D9), A7b Q1b (window courses only class A), A7c Q2 head-to-head."""
    F1 = FR["Q1"]
    rows = O.contrast(F1, "Q1", "A7a", "Q1 by UTI type (D9)", expo="q1_type",
                      levels=["single_lower", "recurrent_lower_only", "recurrent_any_upper", "single_upper"])
    rows += O.contrast(F1.loc[F1["q1b_only_A"]], "Q1", "A7b", "Q1b: window courses only nitrofurantoin/fosfomycin")
    for a, b in (("D", "C"), ("D", "B"), ("B", "C")):
        rows += O.contrast(FR["Q2"], "Q2", "A7c", "Q2 head-to-head", levels=[a, b])
    O.show_rows("secondary", "Secondary contrasts: UTI type, Q1b, Q2 head-to-head (IPTW Cox, MI)", rows,
                note=f"MI m = {K.ANALYSIS['mice_m_sensitivity']} for secondary analyses (O16)")
C.part("5 secondary (A7a-A7c)", p_secondary)
C.end()
##### END F05 #####
