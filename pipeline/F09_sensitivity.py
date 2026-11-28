##### BEGIN F09 sensitivity, bias and validation analyses #####
# Protocol 10 and 12.2 analyses 10-12 and 14-18, 22 (ANALYSIS-MAP A10-A12, A14-A18, A22). Each analysis changes one
# thing from the primary and is fitted with MI m = K.ANALYSIS["mice_m_sensitivity"] (O16). Runs only after the freeze.
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import analysis as A
import cohort as B
import rutirc as R
import outcomes as O

C.begin("F09", "sensitivity, bias and validation analyses")
O.freeze_check()
D = O.load_D()
FR = {q: C.load_state(f"frame_{q}") for q in ("Q1", "Q2")}
QS = ("Q1", "Q2")
DAY = R.DAY


def both(fn, aid, label, **kw):
    rows = []
    for q in QS:
        F = fn(q)
        if F is not None:
            rows += O.contrast(F, q, aid, label, **kw)
    return rows


def p_outcome():
    """A10: outcome confirmed + possible; outcome dated at detection; nephroureterectomy not confirming (2026-11-28)."""
    rows = both(lambda q: B.build(D, B.default_spec(q, outcome_level="confirmed_or_possible"))[0], "A10", "RCC confirmed or possible")
    det = D["rcc"].assign(dx_date=D["rcc"]["detection_date"].fillna(D["rcc"]["dx_date"]))
    rows += both(lambda q: B.build(dict(D, rcc=det), B.default_spec(q))[0], "A10", "RCC dated at detection")
    sens = C.load_state("rcc_without_nephroureterectomy")
    rows += both(lambda q: B.build(dict(D, rcc=sens), B.default_spec(q))[0], "A10", "nephroureterectomy not confirming")
    O.show_rows("outcome", "Outcome definition sensitivity (A10)", rows)
C.part("1 outcome definition (A10)", p_outcome)


def p_uti():
    """A11: episode window 14 days; code-only UTI (no confirmation); culture-confirmed UTI only."""
    codes, rec = C.load_state("uti_codes"), C.load_state("abx_records")
    abx = rec.loc[rec["urinary_active"], ["person_id", "start_date"]].rename(columns={"start_date": "date"}).drop_duplicates()
    cul = pd.concat([D["px"].loc[D["px"]["set_key"].eq("URINE_CULTURE"), ["person_id", "date"]],
                     D["labs"].loc[D["labs"]["set_key"].eq("URINE_CULTURE_LAB"), ["person_id", "date"]]]).drop_duplicates()
    variants = {"episode window 14 days": R.build_episodes(codes, abx, cul, gap_days=14),
                "code-only UTI": R.build_episodes(codes, abx, cul, require_confirmation=False)}
    E = D["episodes"]
    variants["culture-confirmed UTI only"] = E.assign(qualifying=E["confirmed_culture"].astype(bool))
    rows = []
    for lab, Ev in variants.items():
        rows += both(lambda q: B.build(dict(D, episodes=Ev), B.default_spec(q))[0], "A11", lab)
    O.show_rows("uti", "UTI definition sensitivity (A11)", rows)
C.part("2 UTI definition (A11)", p_uti)


def well_observed(F):
    """A12: >= 5 years of EHR history before T0 and >= 1 EHR visit in every complete year of follow-up."""
    v = D["visits"].loc[D["visits"]["is_ehr"], ["person_id", "visit_start_date"]].merge(F[["person_id", "S", "end"]], on="person_id")
    v = v.loc[(v["visit_start_date"] >= v["S"]) & (v["visit_start_date"] < v["end"])]
    v["yr"] = ((v["visit_start_date"] - v["S"]).dt.days // 365.25).astype(int)
    have = v.groupby("person_id")["yr"].nunique()
    full = np.floor(F["t_years"]).astype(int)
    ok = (F["ehr_years_before_t0"] >= K.ANALYSIS["well_observed_history_years"]) & \
         (F["person_id"].map(have).fillna(0).astype(int) >= full)
    return F.loc[ok]


def p_obs():
    """A12 well-observed people; A15 ICD-10 era only (T0 >= 2018-10-01)."""
    rows = both(lambda q: well_observed(FR[q]), "A12", "well observed (>= 5 y history, visit every follow-up year)")
    rows += both(lambda q: FR[q].loc[FR[q]["T0"] >= pd.Timestamp(K.ICD10_ERA_ONLY_T0)], "A15", "ICD-10 era only (T0 >= 2018-10-01)")
    O.show_rows("obs", "Well-observed restriction (A12) and ICD-10 era (A15)", rows)
C.part("3 well observed, ICD-10 era (A12, A15)", p_obs)


def site_complete(F):
    """Persons whose site had EHR visits but no drug records in any calendar year from T0 - 3 y to T0 are left out."""
    sy = C.load_state("site_year", required=False)
    if sy is None:
        return None
    bad = sy.loc[(sy["n_visits"] > 0) & (sy["n_drugs"] == 0), ["src_id", "yr"]]
    m = F[["person_id", "site", "T0"]].merge(bad, left_on="site", right_on="src_id")
    m = m.loc[(m["yr"] >= m["T0"].dt.year - 3) & (m["yr"] <= m["T0"].dt.year)]
    return F.loc[~F["person_id"].isin(m["person_id"])]


def p_site():
    """A14: Cox not stratified by site; leave-one-site-out range; site-year drug completeness restriction."""
    rows = both(lambda q: FR[q], "A14", "not stratified by site", strata=False)
    rows += both(lambda q: site_complete(FR[q]), "A14", "sites-years with drug data only")
    O.show_rows("site", "Site sensitivity (A14)", rows)
    out = []
    for q in QS:
        F = FR[q]
        hr = {lev: [] for lev in O.LEVELS[q][1:]}
        for s in sorted(F["site"].unique()):
            for r in O.contrast(F.loc[F["site"].ne(s)], q, "A14", "leave one site out", method="category"):
                if r.get("status") == "shown":
                    hr[r["level"]].append(r["HR"])
        for lev, v in hr.items():
            out.append({"question": q, "contrast": f"{O.LEVEL_NAMES[lev]} vs {O.LEVEL_NAMES[O.LEVELS[q][0]]}",
                        "sites_left_out": len(v), "HR_min": min(v) if v else np.nan,
                        "HR_median": float(np.median(v)) if v else np.nan, "HR_max": max(v) if v else np.nan})
    C.show("loso", "Leave-one-site-out hazard ratios (A14; missing-category data, IPTW)", pd.DataFrame(out), counts=[],
           note="one fit per EHR site left out; sites are not named; fits with an arm <= 20 events are skipped")
C.part("4 site (A14)", p_site)


def censor_at(F, dates):
    """End follow-up at the first date in `dates` (person -> date) after outcome counting starts; drop people with
    that date between T0 and outcome-counting start (as the cohort builder does for ESRD)."""
    c = F["person_id"].map(dates)
    F = F.loc[~(c < F["S"])].copy()
    c = c.loc[F.index]
    end = F["end"].where(~(c < F["end"]), c)
    F["event"] = (F["dx_date"].notna() & (F["dx_date"] <= end)).astype(int)
    F["death_competing"] = ((F["event"] == 0) & F["death_date"].notna() & (F["death_date"] <= end)).astype(int)
    F["end"], F["t_years"] = end, (end - F["S"]).dt.days / 365.25
    return F.loc[F["t_years"] > 0]


def hematuria_free(F):
    ep = D["episodes"].loc[D["episodes"]["hematuria"].astype(bool), ["person_id", "ep_start"]].merge(F[["person_id", "T0"]], on="person_id")
    h = ep.loc[(ep["ep_start"] < ep["T0"]) & (ep["ep_start"] >= ep["T0"] - 730 * DAY), "person_id"]
    return F.loc[~F["hematuria_2y"] & ~F["person_id"].isin(h)]


def p_misc():
    """A16: excluding kidney stones; excluding hematuria (2 y before T0, incl. cystitis-with-hematuria episodes);
    outpatient index antibiotic only (Q2); censoring at a new non-kidney cancer; Fine-Gray subdistribution HR."""
    rows = both(lambda q: FR[q].loc[~FR[q]["stones"]], "A16", "excluding kidney stones")
    rows += both(lambda q: hematuria_free(FR[q]), "A16", "excluding hematuria before T0")
    rows += O.contrast(FR["Q2"].loc[FR["Q2"]["episode_setting"].eq("outpatient")], "Q2", "A16", "outpatient index antibiotic only")
    oc = D["dx"].loc[D["dx"]["set_key"].eq("OTHER_CANCER")].groupby("person_id")["date"].min()
    rows += both(lambda q: censor_at(FR[q], oc), "A16", "censored at a new non-kidney cancer")
    fg = lambda d: A.fine_gray_rows(d.assign(cause=np.where(d["event"].eq(1), 1, np.where(d["death_competing"].eq(1), 2, 0))))
    rows += both(lambda q: FR[q], "A16", "Fine-Gray subdistribution HR (death competing)", rows_fn=fg)
    O.show_rows("misc", "Other sensitivity analyses (A16)", rows, note="Fine-Gray rows: HR column = subdistribution HR")
C.part("5 other sensitivity (A16)", p_misc)


def p_detect():
    """A17: adding time-varying abdominal imaging after T0; stratified by pyelonephritis (UTI type)."""
    rows = both(lambda q: FR[q], "A17", "adjusted for time-varying imaging after T0", rows_fn=O.imaging_rows_fn(D),
                extra=["imaging_post_t0"])
    F1 = FR["Q1"]
    for t in ("lower_only", "any_upper"):
        rows += O.contrast(F1.loc[F1["uti_type_window"].eq(t)], "Q1", "A17", f"Q1 within {t.replace('_', ' ')} UTI")
    F2 = FR["Q2"]
    for t in ("lower", "upper"):
        rows += O.contrast(F2.loc[F2["episode_type"].eq(t)], "Q2", "A17", f"Q2 within {t} index episode")
    O.show_rows("detect", "Detection-bias sensitivity (A17)", rows)
C.part("6 detection (A17)", p_detect)


def p_negative():
    """A18: negative-control outcomes (D14) with the primary models; negative-control exposure (antifungal and
    antiviral courses) in the Q2b model. Expected HR about 1."""
    rows = []
    for key in K.DECISIONS["D14_nco"]:
        nco = D["dx"].loc[D["dx"]["set_key"].eq(key)].groupby("person_id")["date"].min().rename("dx_date").reset_index()
        nco = nco.assign(level="confirmed", detection_date=pd.NaT)
        rows += both(lambda q: B.build(dict(D, rcc=nco), B.default_spec(q))[0], "A18", f"negative-control outcome {key}")
    nc = C.load_state("nc_courses", required=False)
    if nc is not None:
        rows += [r for r in O.q2b(FR["Q1"], D, K.DECISIONS["D11_lag_years"], "A18", "negative-control exposure in Q2b", nc=nc)
                 if r["contrast"].startswith("nc_cum")]
    O.show_rows("negative", "Negative controls (A18)", rows,
                note="negative-control outcomes: people with the control condition before outcome counting starts are left out")
C.part("7 negative controls (A18)", p_negative)


def p_valid():
    """A22 (descriptive): how confirmed RCC events in the analysed cohorts were confirmed; how Q1 index episodes were
    confirmed. Survey self-report of kidney cancer needs the Personal Medical History items, not extracted yet."""
    rows = []
    rcc = D["rcc"].set_index("person_id")
    for q in QS:
        ev = FR[q].loc[FR[q]["event"].eq(1), "person_id"]
        route = ev.map(rcc["confirm_route"]).value_counts()
        for k in ("procedure", "two_codes", "drug"):
            rows.append({"question": q, "item_group": "RCC events", "item": f"confirmed by {k}", "n": int(route.get(k, 0)),
                         "of_n": len(ev)})
    E1 = FR["Q1"]
    for k in ("confirmed_abx", "confirmed_culture"):
        if k in E1:
            rows.append({"question": "Q1", "item_group": "index UTI episodes", "item": k.replace("_", " by "),
                         "n": int(E1[k].astype(bool).sum()),
                         "of_n": len(E1)})
    T = pd.DataFrame(rows)
    T["share"] = T["n"] / T["of_n"].where(T["of_n"] > 0)
    C.show("valid", "Validation (A22): confirmation routes of RCC events and of index UTI episodes", T,
           counts=["n", "of_n"], deps={"share": ["n", "of_n"]}, group_cols=["question", "item_group"],
           note="confirmation routes are hierarchical (procedure, then two codes, then drug). Survey self-report "
                "concordance: not run (Personal Medical History items not extracted).")
C.part("8 validation (A22)", p_valid)
C.end()
##### END F09 #####
