"""cohort.py - one cohort builder for every analysis (ANALYSIS-STANDARD rule: one build(spec)).

build(D, spec) takes the F02/F03 state frames (dict D) and a spec:
    spec = dict(question="Q1" | "Q2", latency_years=..., outcome_level=..., index_choice=..., ...)
and returns (frame, flow). The frame is participant-level (INTERNAL); flow is [(label, n persons)] for flow_mask.

Order (P§2, §3, §4, §7): candidates -> T0 window -> age -> EHR history -> data availability -> exclusions ->
no kidney cancer before T0 -> observed at T0 -> [Q2: washout] -> pick one index per person -> latency
(outcome counting starts T0 + latency) -> follow-up and outcome -> covariates.
The index is picked BEFORE the latency step, so the choice of a person's index never depends on later outcomes.
"""
import numpy as np
import pandas as pd

import _config as K
import rutirc as R

DAY = R.DAY
STUDY_CLASSES_Q2 = {"A": "A", "B": "B", "C": "C", "D1": "D", "D2": "D", "D3": "D"}   # P§2.2 (D = beta-lactam)


def default_spec(question, **kw):
    s = dict(question=question, latency_years=K.DECISIONS["D11_latency_years"], lag_years=K.DECISIONS["D11_lag_years"],
             outcome_level=K.DECISIONS["D12_outcome_level"], index_choice=K.IMPLEMENTATION["O1_index_choice"],
             censor_esrd=K.DECISIONS["D20_censor_esrd"], hereditary_before_t0_only=K.DECISIONS["D6_hereditary_before_t0_only"],
             t0_window=K.DECISIONS["D5_t0_window"], age_min=K.AGE_MIN, name=question)
    s.update(kw)
    return s


def spec_text(s):
    return (f"{s['name']}: T0 {s['t0_window'][0]}..{s['t0_window'][1]}, age>={s['age_min']}, latency {s['latency_years']} y, "
            f"outcome {s['outcome_level']}, index {s['index_choice']}, ESRD censor {s['censor_esrd']}")


# -----------------------------------------------------------------------------------------------------------
def _first_before(ev, idx, col, keys=None, strict=True):
    """Per candidate row of idx: any event (person_id, date[, set_key in keys]) before idx[col]."""
    e = ev if keys is None else ev.loc[ev["set_key"].isin(keys)]
    m = idx[["person_id", col]].reset_index().merge(e[["person_id", "date"]], on="person_id")
    hit = m.loc[(m["date"] < m[col]) if strict else (m["date"] <= m[col]), "index"].unique()
    return pd.Series(idx.index.isin(hit), index=idx.index)


def _count_in(ev, idx, lo_col, hi_col, date_col="date", distinct=True):
    """Per candidate: number of (distinct) event dates with lo <= date < hi."""
    m = idx[["person_id", lo_col, hi_col]].reset_index().merge(ev[["person_id", date_col]], on="person_id")
    m = m.loc[(m[date_col] >= m[lo_col]) & (m[date_col] < m[hi_col])]
    n = m.groupby("index")[date_col].nunique() if distinct else m.groupby("index").size()
    return idx.index.to_series().map(n).fillna(0).astype(int)


def _exclusions_before(D, idx, col, spec):
    """P§4 exclusions measured in all history before T0. Returns dict label -> boolean (True = excluded)."""
    dx = D["dx"]
    out = {}
    out["prior kidney cancer (codes, history codes or kidney-cancer drugs)"] = (
        _first_before(dx, idx, col, ["RCC", "RENAL_PELVIS_CA", "HX_KIDNEY_CA"]) | _first_before(D["rcc_drug_starts"], idx, col))
    out["other prior cancer except non-melanoma skin"] = _first_before(dx, idx, col, ["OTHER_CANCER"])
    esrd = D["esrd"].rename(columns={"esrd_date": "date"})
    out["ESRD, dialysis or kidney transplant"] = _first_before(esrd, idx, col)
    if spec["hereditary_before_t0_only"]:
        out["VHL, tuberous sclerosis or polycystic kidney disease"] = _first_before(dx, idx, col, ["VHL", "TSC", "PKD"])
    else:
        ever = set(dx.loc[dx["set_key"].isin(["VHL", "TSC", "PKD"]), "person_id"])
        out["VHL, tuberous sclerosis or polycystic kidney disease"] = idx["person_id"].isin(ever)
    rep = R.repeat_codes_apart(dx.loc[dx["set_key"].eq("CATHETER_REPEAT")], 90).rename("date").reset_index()
    out["chronic catheter or urinary diversion"] = _first_before(dx, idx, col, ["CATHETER_ANY"]) | _first_before(rep, idx, col)
    out["neurogenic bladder"] = _first_before(dx, idx, col, ["NEUROGENIC_BLADDER"])
    out["congenital urologic anomaly"] = _first_before(dx, idx, col, ["CONGENITAL_URO"])
    return out


def _data_available(D, idx, col):
    """P§3.3 person level: >= 1 EHR condition record and >= 1 EHR drug record in the 3 years before T0 (month level)."""
    am = D["activity_months"]
    m = idx[["person_id", col]].reset_index().merge(am, on="person_id")
    m = m.loc[(m["mon"] >= m[col] - K.HISTORY_DAYS * DAY - 31 * DAY) & (m["mon"] < m[col])]
    g = m.groupby("index")[["has_cond", "has_drug"]].max()
    ok = (g["has_cond"] > 0) & (g["has_drug"] > 0)
    return idx.index.to_series().map(ok).fillna(False).astype(bool)


# -----------------------------------------------------------------------------------------------------------
def candidates_q1(D):
    E = D["episodes"]
    cand = R.index_candidates(E)
    cand = R.recurrence_in_window(E, cand)
    cand["T0"] = cand["E1"] + 365 * DAY
    cand["ref"] = cand["E1"]                         # P§9.2: Q1 confounders measured before E1
    return cand


def candidates_q2(D):
    """P§2.2: qualifying episodes treated with exactly one study class started within +/-3 d of the first code.
    T0 = the antibiotic course start."""
    E = D["episodes"].loc[D["episodes"]["qualifying"]].reset_index(drop=True)
    cs = D["courses"]
    w = K.DECISIONS["D8_confirm_abx_days"]
    scope = list(STUDY_CLASSES_Q2) if K.IMPLEMENTATION["O4_q2_multi_class_scope"] != "A-G" else list(K.ABX_CLASSES)
    m = E[["person_id", "ep_start", "ep_type"]].reset_index().merge(cs.loc[cs["abx_class"].isin(scope)], on="person_id")
    m = m.loc[(m["c_start"] >= m["ep_start"] - w * DAY) & (m["c_start"] <= m["ep_start"] + w * DAY)]
    g = m.groupby("index").agg(n_classes=("abx_class", "nunique"), abx_class=("abx_class", "first"),
                               T0=("c_start", "min"), setting=("setting", "first"))
    g = g.loc[g["n_classes"].eq(1) & g["abx_class"].isin(list(STUDY_CLASSES_Q2))]
    cand = E.loc[g.index, ["person_id", "ep_start", "ep_type"]].join(g)
    cand["q2_class"] = cand["abx_class"].map(STUDY_CLASSES_Q2)
    cand["ref"] = cand["T0"]
    return cand.rename(columns={"setting": "episode_setting", "ep_type": "episode_type"}).reset_index(drop=True)


def build(D, spec):
    """Returns (frame, flow). frame: one row per analysed person with exposure, time, event and covariates."""
    q = spec["question"]
    cand = candidates_q1(D) if q == "Q1" else candidates_q2(D)
    flow = [("persons with a qualifying index candidate" + (" (clean 365-d window)" if q == "Q1" else
             " (episode treated with one study class)"), cand["person_id"].nunique())]
    keep = pd.Series(True, index=cand.index)

    def step(label, ok):
        nonlocal keep
        keep &= ok.reindex(cand.index).fillna(False).astype(bool)
        flow.append((label, cand.loc[keep, "person_id"].nunique()))

    lo, hi = (pd.Timestamp(x) for x in spec["t0_window"])
    step(f"T0 between {spec['t0_window'][0]} and {spec['t0_window'][1]}", cand["T0"].between(lo, hi))
    bd = R.birth_dates(D["persons"])
    import aou_core as C
    age = C.completed_age(R.dmap(cand["person_id"], bd), cand["T0"])
    cand["age_t0"] = age.values
    step(f"age >= {spec['age_min']} at T0", pd.Series(age.values >= spec["age_min"], index=cand.index))
    visit_days = D["visits"].loc[D["visits"]["is_ehr"], ["person_id", "visit_start_date"]].rename(
        columns={"visit_start_date": "date"}).drop_duplicates()
    step(">= 3 years of EHR history (visit in each of the 3 years before T0)", R.history_ok(cand, visit_days, D["ehr_span"]))
    step("diagnosis and drug data in the 3 years before T0", _data_available(D, cand, "T0"))
    for label, exc in _exclusions_before(D, cand, "T0", spec).items():
        step("no " + label + " before T0", ~exc)
    span = D["ehr_span"].set_index("person_id")
    death = D["death"].set_index("person_id")["death_date"]
    obs = (R.dmap(cand["person_id"], span["last_ehr"]) >= cand["T0"]) & ~(R.dmap(cand["person_id"], death) < cand["T0"])
    step("alive with EHR activity on or after T0", obs)
    if q == "Q2":
        cs = D["courses"]
        m = cand[["person_id", "T0"]].reset_index().merge(cs[["person_id", "c_start"]], on="person_id")
        wash = m.loc[(m["c_start"] < m["T0"]) & (m["c_start"] >= m["T0"] - K.WASHOUT_Q2_DAYS * DAY), "index"].unique()
        step("no study antibiotic in the 90 days before T0 (new user)", pd.Series(~cand.index.isin(wash), index=cand.index))
    # one index per person (O1)
    elig = cand.loc[keep].sort_values(["person_id", "T0"])
    if spec["index_choice"] == "first_ever":
        first = cand.sort_values(["person_id", "T0"]).groupby("person_id").head(1).index
        elig = elig.loc[elig.index.isin(first)]
    F = elig.groupby("person_id").head(1).reset_index(drop=True)
    if q == "Q1":
        flow.append(("in the recurrent or single group (intermediate reported separately)",
                     int(F["q1_group"].isin(["recurrent", "single"]).sum())))
    # outcome, latency, follow-up (P§7, §8)
    O = D["rcc"]
    if spec["outcome_level"] == "confirmed":
        O = O.loc[O["level"].eq("confirmed")]
    dx_date = R.dmap(F["person_id"], O.set_index("person_id")["dx_date"])
    F["S"] = F["T0"] + pd.to_timedelta(365.25 * spec["latency_years"], unit="D").round("D")
    esrd = R.dmap(F["person_id"], D["esrd"].set_index("person_id")["esrd_date"]) if spec["censor_esrd"] else pd.Series(pd.NaT, index=F.index)
    F["death_date"] = R.dmap(F["person_id"], death)
    F["last_ehr"] = R.dmap(F["person_id"], span["last_ehr"])
    cutoff = pd.Timestamp(K.FOLLOWUP_CUTOFF)
    pre = (dx_date < F["S"]) | (F["death_date"] < F["S"]) | (F["last_ehr"] < F["S"]) | (esrd < F["S"]) | (F["S"] >= cutoff)
    base = F["q1_group"].isin(["recurrent", "single", "intermediate"]) if q == "Q1" else pd.Series(True, index=F.index)
    F = F.loc[~pre].copy()
    dx_date, esrd = dx_date.loc[F.index], esrd.loc[F.index]
    flow.append((f"no RCC, death, ESRD or end of EHR before outcome counting starts (T0 + {spec['latency_years']} y)",
                 int(F["q1_group"].isin(["recurrent", "single"]).sum()) if q == "Q1" else len(F)))
    end = R.row_min([dx_date, F["death_date"], F["last_ehr"], esrd, pd.Series(cutoff, index=F.index)])
    F["end"] = end
    F["event"] = (dx_date.notna() & (dx_date <= end)).astype(int)
    F["dx_date"] = dx_date
    F["t_years"] = (F["end"] - F["S"]).dt.days / 365.25
    F["death_competing"] = ((F["event"] == 0) & F["death_date"].notna() & (F["death_date"] <= F["end"])).astype(int)
    enroll = R.dmap(F["person_id"], D["enroll"].set_index("person_id")["enroll_date"])
    F["entry_years"] = ((enroll - F["S"]).dt.days / 365.25).clip(lower=0).fillna(0)
    F["years_t0_to_enrollment"] = (enroll - F["T0"]).dt.days / 365.25
    F = add_covariates(D, F, q)
    F.attrs["spec"] = spec_text(spec)
    return F, flow


# -----------------------------------------------------------------------------------------------------------
def _cat(series, bins, labels):
    return pd.cut(series, bins=bins, labels=labels, right=False).astype(object).where(series.notna(), "Missing")


def add_covariates(D, F, q):
    """P§9.2, §9.4. 'ref' = E1 for Q1 (before the exposure window), T0 for Q2."""
    dx, px = D["dx"], D["px"]
    F = F.copy()
    F["ref_minus_3y"] = F["ref"] - K.HISTORY_DAYS * DAY
    F["ref_minus_1y"] = F["ref"] - 365 * DAY
    F["t0_year"] = F["T0"].dt.year
    F["coding_era"] = np.where(F["T0"] < pd.Timestamp(K.ICD10_START), "ICD-9", "ICD-10")
    # sex: survey sex at birth, else EHR gender (P§9.3)
    sv = D["survey"]
    sx = sv.loc[sv["question_code"].isin([k for k, v in K.SURVEY_ITEMS.items() if v == "sex_at_birth"])]
    sx = sx.assign(level=sx["answer_code"].map(K.CROSSWALKS["sex_at_birth"]))
    sx = sx.loc[sx["level"].isin(["Male", "Female"])].groupby("person_id")["level"].first()
    gender = D["persons"].set_index("person_id")["gender"].astype(str).str.title()
    F["sex"] = F["person_id"].map(sx).fillna(F["person_id"].map(gender).where(lambda s: s.isin(["Male", "Female"]))).fillna("Missing")
    # race / ethnicity (explicit crosswalk, multi-select; Hispanic separate; M9, M10)
    re_ = sv.loc[sv["question_code"].isin([k for k, v in K.SURVEY_ITEMS.items() if v == "race_eth"])]
    re_ = re_.assign(level=re_["answer_code"].map(K.CROSSWALKS["race_eth"]).fillna("Not mapped"))
    sets = re_.groupby("person_id")["level"].agg(lambda s: set(s))

    def race_eth(s):
        if not isinstance(s, set):
            return "Not answered"
        if "Hispanic" in s:
            return "Hispanic"
        real = s - {"Not answered", "Not mapped"}
        if len(real) > 1:
            return "More than one"
        if len(real) == 1:
            return "Non-Hispanic " + next(iter(real))
        return "Not mapped" if "Not mapped" in s else "Not answered"
    F["race_eth"] = F["person_id"].map(sets).map(race_eth)
    F["hispanic"] = F["race_eth"].eq("Hispanic")
    pr = sv.loc[sv["answer_name"].fillna("").str.contains("puerto", case=False, regex=True), "person_id"].unique()
    F["puerto_rican"] = F["person_id"].isin(pr)
    for var in ("education", "income"):
        qq = [k for k, v in K.SURVEY_ITEMS.items() if v == var]
        a = sv.loc[sv["question_code"].isin(qq)]
        cw = K.CROSSWALKS.get(var)
        lev = a["answer_code"].map(cw) if cw else a["answer_name"]
        lev = lev.where(~a["answer_code"].isin(K.SKIP_CODES), "Not answered")
        F[var] = F["person_id"].map(pd.Series(lev.values, index=a["person_id"].values).groupby(level=0).first()).fillna("Not answered")
    # smoking: EHR codes before ref, or survey 100 cigarettes = yes (P§9.4)
    sm_q = [k for k, v in K.SURVEY_ITEMS.items() if v == "smoke_100"]
    sm_yes = sv.loc[sv["question_code"].isin(sm_q) & sv["answer_name"].fillna("").str.contains(r"^yes", case=False, regex=True), "person_id"]
    F["smoking"] = np.where(_first_before(dx, F, "ref", ["SMOKING_DX"]) | F["person_id"].isin(sm_yes), "Ever", "Never/unknown")
    # BMI: closest measurement in the 3 years before ref; tie-break latest date then smallest value (M28)
    lab = D["labs"]
    b = lab.loc[lab["set_key"].eq("BMI"), ["person_id", "date", "value_as_number"]].merge(F[["person_id", "ref", "ref_minus_3y"]], on="person_id")
    lo_, hi_ = K.LABS["BMI"]["range"]
    b = b.loc[(b["date"] < b["ref"]) & (b["date"] >= b["ref_minus_3y"]) & b["value_as_number"].between(lo_, hi_)]
    b = b.sort_values(["person_id", "date", "value_as_number"], ascending=[True, False, True]).groupby("person_id").head(1)
    F["bmi"] = F["person_id"].map(b.set_index("person_id")["value_as_number"])
    F["bmi_cat"] = _cat(F["bmi"], [0, 25, 30, 35, 200], ["<25", "25-29.9", "30-34.9", ">=35"])
    # drugs by ATC group (records per person-date)
    drug, atc = D["drug"], D["atc_drugs"]
    dg = drug.merge(atc[["concept_id", "set_key"]], left_on="drug_concept_id", right_on="concept_id")
    dg = dg.rename(columns={"start_date": "date"})
    ah = dg.loc[dg["set_key"].eq("ANTIHYPERTENSIVE")]
    F["hypertension"] = _first_before(dx, F, "ref", ["HYPERTENSION"]) | (_count_in(ah, F, "ref_minus_3y", "ref") >= 2)
    a1c = lab.loc[lab["set_key"].eq("HBA1C") & (lab["value_as_number"] >= 6.5), ["person_id", "date"]]
    F["diabetes"] = (_first_before(dx, F, "ref", ["DIABETES"]) | _first_before(dg.loc[dg["set_key"].eq("GLUCOSE_LOWERING")], F, "ref")
                     | _first_before(a1c, F, "ref"))
    ckd_lab = D["ckd_lab"].rename(columns={"ckd_lab_date": "date"})
    F["ckd_pre_e1"] = _first_before(dx, F, "ref", ["CKD"]) | _first_before(ckd_lab, F, "ref")
    F["stones"] = _first_before(dx, F, "ref", ["STONES"]) | _first_before(px, F, "ref", ["STONE_PX"])
    imm = dg.loc[dg["set_key"].isin(["SYSTEMIC_STEROID", "IMMUNOSUPPRESSANT"])].copy()
    ds = pd.to_numeric(imm["days_supply"], errors="coerce").fillna(K.DECISIONS["D10_default_course_days"]).clip(lower=1)
    imm["dend"] = imm["date"] + pd.to_timedelta(ds, unit="D")
    mi = F[["person_id", "ref", "ref_minus_1y"]].reset_index().merge(imm[["person_id", "date", "dend"]], on="person_id")
    mi["days"] = ((mi[["dend", "ref"]].min(axis=1) - mi[["date", "ref_minus_1y"]].max(axis=1)).dt.days).clip(lower=0)
    imm_days = F.index.to_series().map(mi.groupby("index")["days"].sum()).fillna(0)
    F["immunosuppression"] = (_first_before(dx, F, "ref", ["HIV", "IMMUNO_DX"]) | (imm_days >= 90))
    # healthcare use and detection (P§9.4)
    vd = D["visits"].loc[D["visits"]["is_ehr"], ["person_id", "visit_start_date"]].rename(columns={"visit_start_date": "date"})
    F["visits_prior_year"] = _count_in(vd, F, "ref_minus_1y", "ref")
    F["imaging_prior_3y"] = _count_in(px.loc[px["set_key"].eq("IMAGING_ABD")], F, "ref_minus_3y", "ref")
    ep = D["episodes"].loc[D["episodes"]["qualifying"], ["person_id", "ep_start"]].rename(columns={"ep_start": "date"})
    hi_col = "ref_minus_1y" if q == "Q1" else "ref"             # Q1: years -3 to -1 before E1 (P§9.2)
    F["prior_uti_episodes"] = _count_in(ep, F, "ref_minus_3y", hi_col)
    cs = D["courses"].rename(columns={"c_start": "date"})
    F["prior_abx_courses"] = _count_in(cs, F, "ref_minus_3y", "ref", distinct=False)
    F["prior_abx_cat"] = _cat(F["prior_abx_courses"].astype(float), [0, 1, 2, 6, 10 ** 6], ["0", "1", "2-5", ">=6"])
    F["prophylaxis_before"] = _first_before(D["prophylaxis"].rename(columns={"p_start": "date"}), F, "ref")
    F["hematuria_2y"] = _count_in(dx.loc[dx["set_key"].eq("HEMATURIA")], F.assign(lo2=F["T0"] - 730 * DAY), "lo2", "T0") > 0
    first_ehr = R.dmap(F["person_id"], D["ehr_span"].set_index("person_id")["first_ehr"])
    F["ehr_years_before_t0"] = (F["T0"] - first_ehr).dt.days / 365.25
    F["resistant_org"] = _first_before(dx, F, "ref", ["RESISTANT_ORG"])
    F["pcn_allergy"] = _first_before(dx, F, "ref", ["PCN_ALLERGY"])
    eg = D["egfr"].merge(F[["person_id", "ref", "ref_minus_1y"]], on="person_id")
    eg = eg.loc[(eg["date"] < eg["ref"]) & (eg["date"] >= eg["ref_minus_1y"])].sort_values(
        ["person_id", "date", "egfr"], ascending=[True, False, True]).groupby("person_id").head(1)
    F["egfr_cat"] = _cat(F["person_id"].map(eg.set_index("person_id")["egfr"]), [0, 30, 60, 90, 10 ** 4], ["<30", "30-59", "60-89", ">=90"])
    if q == "Q1":
        # Q1b (P§1.1): urinary courses in the window were only class A
        w = D["courses"].merge(F[["person_id", "E1"]], on="person_id")
        w = w.loc[w["urinary"] & (w["c_start"] >= w["E1"]) & (w["c_start"] < w["E1"] + 365 * DAY)]
        only_a = w.groupby("person_id")["abx_class"].agg(lambda s: set(s) == {"A"})
        F["q1b_only_A"] = F["person_id"].map(only_a).fillna(False).astype(bool)
    for c in ("ref_minus_3y", "ref_minus_1y"):
        del F[c]
    return F
