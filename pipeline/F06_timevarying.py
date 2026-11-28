##### BEGIN F06 Q2b dose-response, lag and latency, delayed entry #####
# Protocol 2.3 Q2b, 12.2 analyses 9 and 13 (ANALYSIS-MAP A7d, A9, A13, A21). Runs only after the freeze.
import pandas as pd
import aou_core as C
import _config as K
import cohort as B
import outcomes as O

C.begin("F06", "Q2b time-varying courses, lag / latency, delayed entry")
O.freeze_check()
D = O.load_D()
FR = {q: C.load_state(f"frame_{q}") for q in ("Q1", "Q2")}
LAG = K.DECISIONS["D11_lag_years"]


def p_q2b():
    """A7d: cumulative courses by class (lagged 2 y), Q1 eligibles including intermediates."""
    rows = O.q2b(FR["Q1"], D, LAG, "A7d", f"Q2b cumulative courses by class, lag {LAG} y")
    O.show_rows("q2b", "Q2b: time-varying cumulative antibiotic courses by class (Cox, adjusted, MI)", rows,
                note="Adjusted for ADJ_Q1 and lagged cumulative UTI episodes; courses counted from E1; categories 0 / 1 / "
                     "2-5 / >=6 (D10); per-course trend capped at 10. Not IPTW (P§2.3).")
C.part("1 Q2b (A7d)", p_q2b)


def p_lag():
    """A9: latency 1 and 3 years (Q1, Q2); lag and latency together and apart (Q2b)."""
    rows = []
    for lat in (1, 3):
        for q in ("Q1", "Q2"):
            F, _ = B.build(D, B.default_spec(q, latency_years=lat))
            rows += O.contrast(F, q, "A9", f"latency {lat} y")
    for lag, lat in ((1, 2), (3, 2), (2, 1), (2, 3), (1, 1), (3, 3)):
        F, _ = B.build(D, B.default_spec("Q1", latency_years=lat))
        rows += O.q2b(F, D, lag, "A9", f"Q2b lag {lag} y, latency {lat} y")
    O.show_rows("lag", "Lag and latency sensitivity (A9)", rows)
C.part("2 lag / latency (A9)", p_lag)


def p_entry():
    """A13: delayed entry at enrollment (left truncation): follow-up counted only from enrollment."""
    rows = []
    for q in ("Q1", "Q2"):
        rows += O.contrast(FR[q], q, "A13", "delayed entry at enrollment", entry="entry_years")
    O.show_rows("entry", "Delayed entry at enrollment (A13)", rows,
                note="persons whose follow-up ends before enrollment are left out of this analysis")
C.part("3 delayed entry (A13)", p_entry)


def p_q2c():
    """A21: Q2c needs a cohort of people with no UTI ever, which F02 does not extract (candidates = UTI or RCC codes)."""
    C.show("q2c", "Q2c (A21): not run", pd.DataFrame([{"analysis_id": "A21", "status":
           "not run: needs a no-UTI extraction (F02 candidates are people with a UTI or RCC code)"}]), counts=[])
C.part("4 Q2c (A21)", p_q2c)
C.end()
##### END F06 #####
