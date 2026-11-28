"""outcomes.py - the RUTI-RCC outcome-model workflow shared by F05-F10 (one engine, like cohort.build).

Rules built in here, so no cell can skip them:
  * freeze_check(): every outcome cell stops unless the blind F00-F04 output is saved in workbench_outputs/,
    named in _config.RUN["blind_frozen_output"], and its F04 code hash is the F04 that built this run's cohorts.
  * an HR is shown only when each compared arm has >= 21 events (D15 per-cell rule, applied everywhere); event
    counts are never printed outside the F05 rates table (so nested analyses cannot be differenced).
  * missing data: "mice" (primary, P§12.2 #8), "complete_case" or "category" (missing as its own level).
Participant-level frames stay inside the Workbench; everything returned is aggregate.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

import aou_core as C
import _config as K
import analysis as A
import cohort as B

MIN_EV = K.ANALYSIS["min_events_per_arm"]
EXPO = {"Q1": "q1_group", "Q2": "q2_class"}
LEVELS = {"Q1": ["single", "recurrent"], "Q2": ["A", "B", "C", "D"]}
LEVEL_NAMES = {"single": "single UTI", "recurrent": "recurrent UTI", "intermediate": "intermediate",
               "A": "nitrofurantoin/fosfomycin", "B": "TMP-SMX/trimethoprim", "C": "fluoroquinolone", "D": "beta-lactam",
               "single_lower": "single lower UTI", "single_upper": "single upper UTI",
               "recurrent_lower_only": "recurrent, lower only", "recurrent_any_upper": "recurrent, any upper"}
NOT_SHOWN = "not shown: an arm has <=20 events"


def adj(q):
    return list(K.ADJ_Q1 if q == "Q1" else K.ADJ_Q2)


# =============================================================================
# FREEZE GATE
# =============================================================================
def repo_root():
    return Path(__file__).resolve().parent.parent


def stamp_f04(code_hash):
    """Called at the end of F04: records which code built this run's cohorts."""
    p = C._run() / "state" / "F04_stamp.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"code": code_hash, "version": C.CFG.get("version", "")}))


def freeze_check():
    """STUDY-LIFECYCLE stage 8 / P§12.1: outcome models run only after the blind cohort output is frozen."""
    name = (K.RUN.get("blind_frozen_output") or "").strip()
    stop = "STOPPED (freeze gate, expected before the blind output is frozen): "
    if not name:
        raise RuntimeError(stop + "_config.RUN['blind_frozen_output'] is empty. Run F00-F04, copy out ALL_TABLES.txt, "
                           "save it in workbench_outputs/, review it, then set its file name there and pull.")
    p = Path(name) if Path(name).is_absolute() else repo_root() / "workbench_outputs" / name
    if not p.exists():
        raise RuntimeError(stop + f"{p} not found (git pull?).")
    txt = p.read_text()
    if "--- F04.events_total" not in txt or "--- F04.gate" not in txt:
        raise RuntimeError(stop + f"{p.name} has no F04 events_total / gate tables.")
    hashes = set(re.findall(r"code ([0-9a-f]{12})", txt[txt.index("--- F04."):]))
    sp = C._run() / "state" / "F04_stamp.json"
    if not sp.exists():
        raise RuntimeError(stop + f"run folder {C._run().name} has no F04 stamp: F05+ must run in the run folder that "
                           "produced the frozen output.")
    h = json.loads(sp.read_text())["code"]
    if h not in hashes:
        raise RuntimeError(stop + f"this run's F04 used code {h}; the frozen output has {sorted(hashes)}.")
    print(f"freeze gate passed: {p.name}, F04 code {h}", flush=True)


# =============================================================================
# DATA FOR MODELS
# =============================================================================
def load_D():
    names = ["dx", "px", "drug", "visits", "labs", "persons", "death", "enroll", "survey", "ehr_span", "activity_months",
             "episodes", "courses", "prophylaxis", "egfr", "esrd", "ckd_lab", "rcc"]
    D = {n: C.load_state(n) for n in names}
    D["atc_drugs"] = C.load_state("atc_drugs", internal=False)
    G = C.load_state("drug_glossary", internal=False)
    D["rcc_drug_starts"] = D["drug"].merge(G.loc[G["rcc_drug"], ["concept_id"]], left_on="drug_concept_id",
                                           right_on="concept_id")[["person_id", "start_date"]].rename(columns={"start_date": "date"})
    return D


def prep(F):
    """Model-ready copy: log(1 + x) for count covariates."""
    G = F.copy()
    for c in K.ANALYSIS["count_covariates_log1p"]:
        if c in G:
            G[c] = np.log1p(pd.to_numeric(G[c], errors="coerce").astype(float))
    return G


def targets(q):
    return ["bmi", "education", "income"] + (["egfr"] if q == "Q2" else [])


def is_missing(F, c):
    return F[c].isna() if c in ("bmi", "egfr") else F[c].eq("Not answered") | F[c].isna()


def complete_case(F, q):
    ok = pd.Series(True, index=F.index)
    for c in targets(q):
        ok &= ~is_missing(F, c)
    return F.loc[ok]


def impute(F, q, expo, m, seed=20261128):
    """MICE with predictive mean matching (P§12.2 #8): targets bmi, egfr (Q2), education, income (ordinal);
    predictors = exposure, the other adjustment covariates, the event indicator and the Nelson-Aalen cumulative
    hazard (White & Royston 2009). Returns m completed frames."""
    tg = [c for c in targets(q) if c in F]
    T = pd.DataFrame({c: (pd.to_numeric(F[c], errors="coerce") if c in ("bmi", "egfr") else
                          F[c].map({lev: i for i, lev in enumerate(K.ORDINAL_LEVELS[c])})).astype(float).values for c in tg})
    preds = [c for c in adj(q) if c not in tg + ["bmi_cat", "egfr_cat"]]
    G = prep(F).reset_index(drop=True).assign(_expo=F[expo].astype(str).values)
    Xp, _ = C.design(G, preds + ["_expo"])
    Xp = A.drop_constant(Xp)[0]
    Xp["_event"] = F["event"].values
    Xp["_H"] = A.nelson_aalen(F["t_years"].values, F["event"].values)
    Z = A.mice_pmm(T, Xp.values, m=m, iters=K.ANALYSIS["mice_iters"], seed=seed, fit_rows=K.ANALYSIS["mice_fit_rows"])
    out = []
    for z in Z:
        D = F.copy()
        for c in tg:
            if c in ("bmi", "egfr"):
                D[c] = z[c].values
            else:
                lv = np.array(K.ORDINAL_LEVELS[c], dtype=object)
                D[c] = lv[np.clip(np.round(z[c].values).astype(int), 0, len(lv) - 1)]
        out.append(B.recode_imputed(D))
    return out


def datasets(F, q, expo, method, m):
    if method == "mice":
        return impute(F, q, expo, m)
    if method == "complete_case":
        return [complete_case(F, q)]
    if method == "category":
        return [F]
    raise ValueError(method)


def arm_events(F, expo, levels, event="event"):
    return {lev: int(F.loc[F[expo].eq(lev), event].sum()) for lev in levels}


def person_weights(d, expo, levels, covs):
    """Stabilized, truncated IPTW per person (index aligned with d) and PS diagnostics."""
    G = prep(d).reset_index(drop=True)
    X, _ = C.design(G, [c for c in covs if c in G])
    sw, info = A.ps_weights(G[expo], X, levels, K.ANALYSIS["ps_truncate_pct"])
    info["X"] = X
    return pd.Series(sw, index=d.index), info


def expo_matrix(R, expo, levels):
    return pd.DataFrame({f"{expo}[{lev}]": R[expo].eq(lev).astype(float).values for lev in levels[1:]}, index=R.index)


def cox_fit(R, Xcols, w=None, strata=True, entry=None, time="t_years", event="event", cluster=None):
    X = pd.concat(Xcols, axis=1) if isinstance(Xcols, list) else Xcols
    return A.cox(R[time].values, R[event].values, X.values, w=None if w is None else np.asarray(w),
                 entry=None if entry is None else R[entry].values,
                 strata=R["site"].values if strata else None,
                 cluster=None if cluster is None else R[cluster].values, names=list(X.columns))


# =============================================================================
# CONTRASTS
# =============================================================================
def contrast(F, q, aid, label, expo=None, levels=None, method="mice", m=None, covs=None, strata=None, entry=None,
             rows_fn=None, extra=(), time="t_years", event="event"):
    """IPTW cause-specific Cox for levels[1:] vs levels[0]. rows_fn(d) -> counting-process rows (start/stop) with
    the `extra` time-varying columns added to the outcome model. Returns one row per compared level."""
    expo = expo or EXPO[q]
    levels = levels or LEVELS[q]
    covs = adj(q) if covs is None else covs
    strata = K.ANALYSIS["site_stratum"] if strata is None else strata
    m = m or (K.ANALYSIS["mice_m"] if aid in ("A5", "A5b") else K.ANALYSIS["mice_m_sensitivity"])
    F = F.loc[F[expo].isin(levels)].reset_index(drop=True)
    if entry:
        F = F.loc[F[entry] < F[time]].reset_index(drop=True)
    if method == "complete_case":
        F = complete_case(F, q).reset_index(drop=True)
    ev = arm_events(F, expo, levels, event)
    base = {"analysis_id": aid, "analysis": label, "question": q, "missing": method,
            "m": m if method == "mice" else 1}
    rows = [dict(base, contrast=f"{LEVEL_NAMES.get(lev, lev)} vs {LEVEL_NAMES.get(levels[0], levels[0])}", level=lev)
            for lev in levels[1:]]
    ok = {lev: ev[levels[0]] >= MIN_EV and ev[lev] >= MIN_EV for lev in levels[1:]}
    if not any(ok.values()):
        for r in rows:
            r["status"] = NOT_SHOWN
        return rows
    sets = datasets(F, q, expo, "mice" if method == "mice" else "category", m)
    betas, ses, smd_max, ess = {lev: [] for lev in levels[1:]}, {lev: [] for lev in levels[1:]}, {}, {}
    status_bad = set()
    ps_conv = True
    for i, d in enumerate(sets):
        d = d.reset_index(drop=True)
        w, info = person_weights(d, expo, levels, covs)
        ps_conv = ps_conv and info["converged"] if i else info["converged"]
        if i == 0:
            for lev in levels[1:]:
                bt = A.balance(info["X"], d[expo].values, lev, levels[0], w.values)
                small = bt["binary"] & (bt["nz_1"].between(1, 20) | bt["nz_0"].between(1, 20))
                smd_max[lev] = float(bt.loc[~small, "smd_weighted"].abs().max()) if (~small).any() else np.nan
                ess[lev] = (info["ess"][lev], info["ess"][levels[0]])
        if rows_fn is None:
            R, ent, wr, cl = d, entry, w.values, None
        else:
            R = rows_fn(d)
            wr = R["person_id"].map(pd.Series(w.values, index=d["person_id"].values)).values
            if "fg_w" in R:
                wr = wr * R["fg_w"].values
            ent, cl = "start", "person_id"
        Xe = expo_matrix(R, expo, levels)
        Xx = [R[[c for c in extra]].astype(float)] if extra else []
        fit = cox_fit(R, [Xe] + Xx, w=wr, strata=strata, entry=ent, time="stop" if rows_fn else time, event=event,
                      cluster=cl)
        if fit["divergent"] or not fit["converged"]:
            status_bad.add(i)
        for j, lev in enumerate(levels[1:]):
            betas[lev].append(fit["beta"][j]); ses[lev].append(fit["se"][j])
    crude = A.cox(F[time].values, F[event].values, expo_matrix(F, expo, levels).values,
                  entry=None if entry is None else F[entry].values)
    for j, (r, lev) in enumerate(zip(rows, levels[1:])):
        if not ok[lev]:
            r["status"] = NOT_SHOWN
            continue
        if status_bad:
            r["status"] = f"not estimable: model did not converge in {len(status_bad)} of {len(sets)} data sets"
            continue
        b, s, dfree = A.rubin_log(betas[lev], ses[lev])
        r["HR"], r["lo"], r["hi"], r["p"] = A.ratio_ci(b, s, dfree)
        r["beta"], r["se"] = b, s
        r["HR_crude"], r["lo_crude"], r["hi_crude"], _ = A.ratio_ci(crude["beta"][j], crude["se"][j])
        r["max_abs_smd_weighted"] = smd_max.get(lev, np.nan)
        r["ess_exposed"], r["ess_reference"] = (round(x, -1) for x in ess.get(lev, (np.nan, np.nan)))   # nearest 10
        r["ps_converged"] = "yes" if ps_conv else "no"
        r["E_value"] = C.evalue_rr(r["HR"])
        r["E_value_ci"] = C.evalue_ci(r["HR"], r["lo"], r["hi"])
        r["status"] = "shown"
    return rows


def show_rows(tid, title, rows, note=None):
    cols = ["analysis_id", "analysis", "question", "contrast", "HR", "lo", "hi", "p", "HR_crude", "lo_crude",
            "hi_crude", "max_abs_smd_weighted", "ess_exposed", "ess_reference", "ps_converged", "E_value", "E_value_ci",
            "missing", "m", "status"]
    T = pd.DataFrame(rows)
    T = T[[c for c in cols if c in T.columns]]
    record(rows)
    return C.show(tid, title, T, counts=[], max_rows=200,
                  note=(note + " | " if note else "") + "HR = IPTW cause-specific hazard ratio (robust SE), site as "
                  "stratum; crude = unweighted; an HR is shown only if each arm has >= 21 events; event counts by arm "
                  "are in F05.rates only.")


def record(rows):
    """Aggregate results registry (no counts) for F10."""
    new = pd.DataFrame(rows).drop(columns=[c for c in ("level",) if c in pd.DataFrame(rows).columns])
    old = C.load_state("results", internal=False, required=False)
    if old is not None and len(old):
        key = ["analysis_id", "analysis", "question", "contrast"]
        k_old = old[key].astype(str).agg("|".join, axis=1)
        k_new = set(new[key].astype(str).agg("|".join, axis=1))
        new = pd.concat([old.loc[~k_old.isin(k_new)], new], ignore_index=True)
    C.save_state(new, "results", internal=False)


# =============================================================================
# RATES, CUMULATIVE INCIDENCE, PH
# =============================================================================
def rates(F, q, expo=None, levels=None):
    expo = expo or EXPO[q]
    levels = levels or sorted(F[expo].dropna().unique())
    rows = []
    for lev in levels:
        g = F.loc[F[expo].eq(lev)]
        py = float(g["t_years"].sum())
        rows.append({"question": q, "group": LEVEL_NAMES.get(lev, lev), "n": len(g), "events": int(g["event"].sum()),
                     "person_years": round(py, 0), "rate_per_1000py": 1000 * g["event"].sum() / py if py > 0 else np.nan})
    return pd.DataFrame(rows)


def cif(F, q, method="mice", m=None, expo=None, levels=None, covs=None):
    """IPTW Aalen-Johansen CIF (death competing) at K.ANALYSIS['cif_times'] and the difference vs the reference;
    pooled over imputations by Rubin's rules with bootstrap within-variance. Unweighted risk-set and event counts
    decide masking (each arm >= 21 at risk and >= 21 events by that time)."""
    expo = expo or EXPO[q]
    levels = levels or LEVELS[q]
    covs = adj(q) if covs is None else covs
    m = m or K.ANALYSIS["mice_m"]
    at = K.ANALYSIS["cif_times"]
    F = F.loc[F[expo].isin(levels)].reset_index(drop=True)
    if method == "complete_case":
        F = complete_case(F, q).reset_index(drop=True)
    cause = np.where(F["event"].eq(1), 1, np.where(F["death_competing"].eq(1), 2, 0))
    gate = {}
    for lev in levels:
        g = F[expo].eq(lev).values
        aj = A.aalen_johansen(F["t_years"].values[g], cause[g], at=at)
        gate[lev] = ((aj["n_at_risk"] >= MIN_EV) & (aj["events_by"] >= MIN_EV)).values
    if not any(gate[lev].any() for lev in levels):
        return pd.DataFrame([{"question": q, "group": LEVEL_NAMES.get(lev, lev), "years_after_counting_start": x,
                              "status": "not shown: <=20 events or <=20 at risk"} for lev in levels for x in at])
    sets = datasets(F, q, expo, "category" if method == "complete_case" else method, m)
    est = {lev: [] for lev in levels}
    var = {lev: [] for lev in levels}
    dest = {lev: [] for lev in levels[1:]}
    dvar = {lev: [] for lev in levels[1:]}
    B_ = K.ANALYSIS["cif_bootstrap"]
    for i, d in enumerate(sets):
        d = d.reset_index(drop=True)
        w, _ = person_weights(d, expo, levels, covs)
        draws, pts = {}, {}
        cz = np.where(d["event"].eq(1), 1, np.where(d["death_competing"].eq(1), 2, 0))
        for lev in levels:
            g = d[expo].eq(lev).values
            c_ = cz[g]
            pts[lev] = A.aalen_johansen(d["t_years"].values[g], c_, w.values[g], at=at)["cif"].values
            draws[lev] = A.cif_bootstrap(d["t_years"].values[g], c_, w.values[g], at=at, B=B_, seed=1000 * i + 7)
            est[lev].append(pts[lev]); var[lev].append(draws[lev].var(axis=0, ddof=1))
        for lev in levels[1:]:
            dest[lev].append(pts[lev] - pts[levels[0]])
            dvar[lev].append((draws[lev] - draws[levels[0]]).var(axis=0, ddof=1))

    def pool(e, v):
        e, v = np.asarray(e), np.asarray(v)
        mm = len(e)
        tot = v.mean(axis=0) + ((1 + 1 / mm) * e.var(axis=0, ddof=1) if mm > 1 else 0)
        return e.mean(axis=0), np.sqrt(tot)
    rows = []
    for lev in levels:
        pe, ps = pool(est[lev], var[lev])
        for k, x in enumerate(at):
            r = {"question": q, "group": LEVEL_NAMES.get(lev, lev), "years_after_counting_start": x}
            if gate[lev][k]:
                r.update(CIF_pct=100 * pe[k], lo=100 * max(pe[k] - 1.96 * ps[k], 0), hi=100 * (pe[k] + 1.96 * ps[k]),
                         status="shown")
            else:
                r["status"] = "not shown: <=20 events or <=20 at risk"
            if lev != levels[0]:
                de, ds = pool(dest[lev], dvar[lev])
                if gate[lev][k] and gate[levels[0]][k]:
                    r.update(RD_pct_vs_ref=100 * de[k], RD_lo=100 * (de[k] - 1.96 * ds[k]), RD_hi=100 * (de[k] + 1.96 * ds[k]))
            rows.append(r)
    return pd.DataFrame(rows)


def ph_check(F, q, expo=None, levels=None, covs=None):
    """P§12.2 #6 on the missing-category data: Schoenfeld-type score test for each exposure term (beta(t) linear in
    log t), and period-specific HRs 2-5 and > 5 years after T0 with a Wald test of equality."""
    expo = expo or EXPO[q]
    levels = levels or LEVELS[q]
    covs = adj(q) if covs is None else covs
    F = F.loc[F[expo].isin(levels)].reset_index(drop=True)
    w, _ = person_weights(F, expo, levels, covs)
    fit = cox_fit(F, expo_matrix(F, expo, levels), w=w.values, strata=K.ANALYSIS["site_stratum"])
    fit_ph = A.cox(F["t_years"].values, F["event"].values, expo_matrix(F, expo, levels).values, w=w.values,
                   strata=F["site"].values if K.ANALYSIS["site_stratum"] else None, ph=True,
                   names=list(expo_matrix(F, expo, levels).columns))
    cut = K.ANALYSIS["ph_split_years"]
    p1 = F.assign(start=0.0, stop=np.minimum(F["t_years"], cut), ev=np.where(F["t_years"] <= cut, F["event"], 0), period=1)
    p2 = F.loc[F["t_years"] > cut].assign(start=cut, stop=lambda x: x["t_years"], ev=lambda x: x["event"], period=2)
    R = pd.concat([p1, p2], ignore_index=True)
    R = R.loc[R["stop"] > R["start"]].reset_index(drop=True)
    X = pd.concat([expo_matrix(R, expo, levels).mul(R["period"].eq(k).astype(float), axis=0).add_suffix(f"_p{k}")
                   for k in (1, 2)], axis=1)
    fp = cox_fit(R, X, w=R["person_id"].map(pd.Series(w.values, index=F["person_id"].values)).values,
                 strata=K.ANALYSIS["site_stratum"], entry="start", time="stop", event="ev", cluster="person_id")
    rows = []
    L = len(levels) - 1
    for j, lev in enumerate(levels[1:]):
        e_arm = {k: int(R.loc[R["period"].eq(k) & R[expo].eq(lev), "ev"].sum()) for k in (1, 2)}
        e_ref = {k: int(R.loc[R["period"].eq(k) & R[expo].eq(levels[0]), "ev"].sum()) for k in (1, 2)}
        ok = {k: e_arm[k] >= MIN_EV and e_ref[k] >= MIN_EV for k in (1, 2)}
        r = {"question": q, "contrast": f"{LEVEL_NAMES.get(lev, lev)} vs {LEVEL_NAMES.get(levels[0], levels[0])}"}
        if arm_events(F, expo, [levels[0], lev])[lev] >= MIN_EV and arm_events(F, expo, [levels[0]])[levels[0]] >= MIN_EV:
            r["p_schoenfeld_logt"] = float(fit_ph["ph"]["p"].iloc[j])
            r["HR_overall"] = float(np.exp(fit["beta"][j]))
        for k, lab in ((1, "2to5y"), (2, "over5y")):
            if ok[k]:
                r[f"HR_{lab}"] = float(np.exp(fp["beta"][j + (k - 1) * L]))
                r[f"lo_{lab}"] = float(np.exp(fp["beta"][j + (k - 1) * L] - 1.96 * fp["se"][j + (k - 1) * L]))
                r[f"hi_{lab}"] = float(np.exp(fp["beta"][j + (k - 1) * L] + 1.96 * fp["se"][j + (k - 1) * L]))
        if ok[1] and ok[2]:
            r["p_equal_periods"] = A.wald_equal(fp, j, j + L)[1]
        r["status"] = "shown" if "HR_overall" in r else NOT_SHOWN
        rows.append(r)
    return pd.DataFrame(rows)


# =============================================================================
# TIME-VARYING ROWS (Q2b, post-T0 imaging, negative-control courses)
# =============================================================================
def _years_from_S(dates, S):
    return (pd.to_datetime(dates) - pd.to_datetime(S)).dt.days / 365.25


def course_changes(D, F, lag_years, nc=None):
    """Cumulative-count increments for Q2b: study-class courses (A, B, C, D) and other classes started on or after
    E1, qualifying UTI episodes on or after E1, and (nc) negative-control courses, each lagged by lag_years."""
    cs = D["courses"].merge(F[["person_id", "E1", "S"]], on="person_id")
    cs = cs.loc[cs["c_start"] >= cs["E1"]]
    cls = cs["abx_class"].map(lambda c: {"A": "A", "B": "B", "C": "C", "D1": "D", "D2": "D", "D3": "D"}.get(c, "other"))
    lag = pd.to_timedelta(round(365.25 * lag_years), unit="D")
    ch = [pd.DataFrame({"person_id": cs["person_id"].values, "var": ("abx_" + cls).values,
                        "u": _years_from_S(cs["c_start"] + lag, cs["S"]).values})]
    ep = D["episodes"].loc[D["episodes"]["qualifying"]].merge(F[["person_id", "E1", "S"]], on="person_id")
    ep = ep.loc[ep["ep_start"] >= ep["E1"]]
    ch.append(pd.DataFrame({"person_id": ep["person_id"].values, "var": "uti_cum",
                            "u": _years_from_S(ep["ep_start"] + lag, ep["S"]).values}))
    if nc is not None:
        n_ = nc.merge(F[["person_id", "E1", "S"]], on="person_id")
        n_ = n_.loc[n_["c_start"] >= n_["E1"]]
        ch.append(pd.DataFrame({"person_id": n_["person_id"].values, "var": "nc_cum",
                                "u": _years_from_S(n_["c_start"] + lag, n_["S"]).values}))
    return pd.concat(ch, ignore_index=True)


def binned(R, var, bins, labels=None):
    """Dummy columns for cumulative-count categories; bins = lower edges, first = reference."""
    v = R[var].values
    out = {}
    for k in range(1, len(bins)):
        hi = bins[k + 1] if k + 1 < len(bins) else np.inf
        lab = (labels[k] if labels else (f"{bins[k]}" if hi == bins[k] + 1 else
                                         f">={bins[k]}" if not np.isfinite(hi) else f"{bins[k]}-{hi - 1}"))
        out[f"{var}[{lab}]"] = ((v >= bins[k]) & (v < hi)).astype(float)
    return pd.DataFrame(out, index=R.index)


def merge_sparse_levels(F, covs, event="event"):
    """O17 (outcome-model adjustment only): a level of a categorical covariate with no events cannot be estimated
    and makes the whole Cox fit diverge (M12); such levels are merged into the most frequent level (a yes/no
    covariate becomes constant and leaves the model)."""
    F = F.copy()
    for c in covs:
        if c not in F:
            continue
        binary = pd.api.types.is_bool_dtype(F[c]) or (pd.api.types.is_numeric_dtype(F[c]) and F[c].nunique() <= 2)
        if pd.api.types.is_numeric_dtype(F[c]) and not binary:
            continue
        v = F[c].astype(str)
        ev = F.groupby(v)[event].sum()
        empty = ev.index[ev.eq(0)]
        if len(empty) and len(ev) > 1:
            top = v.value_counts().index[0]
            if binary:      # a yes/no covariate with no events in one level: constant, dropped from the design
                F[c] = F[c].where(False, F[c][v.eq(top)].iloc[0])
            else:
                F[c] = v.where(~v.isin(empty), top)
    return F


def trim_bins(R, var, bins, event="event"):
    """Top count categories with no events are merged downward (e.g. 0/1/2-5/>=6 -> 0/1/>=2)."""
    b = list(bins)
    while len(b) > 2 and R.loc[R[var] >= b[-1], event].sum() == 0:
        b = b[:-1]
    return b


def q2b(F, D, lag_years, aid, label, method="mice", m=None, nc=None):
    """P§2.3 Q2b: Cox with time-varying cumulative courses by class (0 / 1 / 2-5 / >= 6, lagged), adjusted for
    baseline confounders (ADJ_Q1) and lagged cumulative UTI episodes; also a per-course trend model.
    nc: negative-control courses (A18) added as one more time-varying count."""
    m = m or K.ANALYSIS["mice_m_sensitivity"]
    F = F.loc[F["q1_group"].isin(["recurrent", "single", "intermediate"])].reset_index(drop=True)
    ch = course_changes(D, F, lag_years, nc)
    classes = ["abx_A", "abx_B", "abx_C", "abx_D"]
    vars_ = classes + ["abx_other", "uti_cum"] + (["nc_cum"] if nc is not None else [])
    sets = datasets(F, "Q1", "q1_group", method, m)
    bins = K.ANALYSIS["q2b_course_bins"]
    res, res_tr, names = [], [], None
    gate = {}
    for i, d in enumerate(sets):
        d = merge_sparse_levels(d.reset_index(drop=True), K.ADJ_Q1)
        R = A.split_cumulative(d, ch, vars_)
        if i == 0:
            vb = {v: trim_bins(R, v, bins) for v in classes + ["abx_other"] + (["nc_cum"] if nc is not None else [])}
            vb["uti_cum"] = trim_bins(R, "uti_cum", [1, 2, 3, 6])
        Xtv = pd.concat([binned(R, v, b) for v, b in vb.items()], axis=1)
        Xb, _ = C.design(prep(R), [c for c in K.ADJ_Q1 if c in R])
        X = pd.concat([Xtv, Xb], axis=1)
        X, _ = A.drop_constant(X)
        X, _ = A.drop_collinear(X)
        if i == 0:
            for c in Xtv.columns:
                gate[c] = int(R.loc[Xtv[c].eq(1), "event"].sum()) if c in X.columns else 0
            ref_ev = {v: int(R.loc[R[v].eq(0) if v != "uti_cum" else R[v].le(1), "event"].sum()) for v in vars_}
        fit = cox_fit(R, X, strata=K.ANALYSIS["site_stratum"], entry="start", time="stop", cluster="person_id")
        names = fit["names"]
        res.append((fit["beta"], fit["se"], fit["divergent"] or not fit["converged"]))
        Xl = pd.concat([R[classes].clip(upper=10).astype(float).add_suffix("_per_course"), Xb,
                        binned(R, "uti_cum", vb["uti_cum"])], axis=1)
        Xl, _ = A.drop_constant(Xl)
        Xl, _ = A.drop_collinear(Xl)
        fl = cox_fit(R, Xl, strata=K.ANALYSIS["site_stratum"], entry="start", time="stop", cluster="person_id")
        res_tr.append((fl["beta"], fl["se"], fl["names"], fl["divergent"] or not fl["converged"]))
    rows = []
    bad = any(r[2] for r in res)
    for c in [c for c in gate if c in names]:
        var = c.split("[")[0]
        j = names.index(c)
        r = {"analysis_id": aid, "analysis": label, "question": "Q2b", "contrast": f"{c} vs {'1' if var == 'uti_cum' else '0'}",
             "missing": method, "m": len(sets)}
        if gate[c] >= MIN_EV and ref_ev[var] >= MIN_EV and not bad:
            b, s, dfree = A.rubin_log([x[0][j] for x in res], [x[1][j] for x in res])
            r["HR"], r["lo"], r["hi"], r["p"] = A.ratio_ci(b, s, dfree)
            r["E_value"], r["E_value_ci"] = C.evalue_rr(r["HR"]), C.evalue_ci(r["HR"], r["lo"], r["hi"])
            r["status"] = "shown"
        else:
            r["status"] = "not estimable: model did not converge" if bad else "not shown: <=20 events in the category or reference"
        rows.append(r)
    for v in classes:
        nm = f"{v}_per_course"
        r = {"analysis_id": aid, "analysis": label + " (per course, capped at 10)", "question": "Q2b",
             "contrast": f"{v}: per additional course", "missing": method, "m": len(sets)}
        js = [x[2].index(nm) if nm in x[2] else None for x in res_tr]
        any_ex = any(gate.get(c, 0) >= MIN_EV for c in gate if c.startswith(v + "["))
        if None not in js and any_ex and ref_ev[v] >= MIN_EV and not any(x[3] for x in res_tr):
            b, s, dfree = A.rubin_log([x[0][j] for x, j in zip(res_tr, js)], [x[1][j] for x, j in zip(res_tr, js)])
            r["HR"], r["lo"], r["hi"], r["p"] = A.ratio_ci(b, s, dfree)
            r["status"] = "shown"
        else:
            r["status"] = "not shown: <=20 events with courses of this class, or not estimable"
        rows.append(r)
    return rows


def imaging_rows_fn(D):
    """A17: time-varying 'any abdominal imaging since T0' (not lagged; detection path)."""
    px = D["px"].loc[D["px"]["set_key"].eq("IMAGING_ABD"), ["person_id", "date"]]

    def fn(d):
        im = px.merge(d[["person_id", "T0", "S"]], on="person_id")
        im = im.loc[im["date"] >= im["T0"]]
        ch = pd.DataFrame({"person_id": im["person_id"].values, "var": "imaging_post_t0",
                           "u": _years_from_S(im["date"], im["S"]).values})
        return A.split_cumulative(d, ch, ["imaging_post_t0"], caps={"imaging_post_t0": 1})
    return fn
