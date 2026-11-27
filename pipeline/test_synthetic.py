"""test_synthetic.py - offline test of the RUTI-RCC pipeline. No All of Us data.

    pip install pandas numpy statsmodels duckdb sqlglot pyarrow
    python test_synthetic.py

Part 1: algorithm unit tests on hand-built cases (each with the answer worked out from the protocol text).
Part 2: a synthetic OMOP CDR in DuckDB; the real BigQuery SQL is transpiled by sqlglot and executed; cells F00-F04
run end to end through the same exec loop as RUN_IN_WORKBENCH.py; then every saved table is scanned for a count
1-20 that escaped suppression. Synthetic data exist only here (M27: no fallback in analysis code).
"""
import hashlib
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aou_core as C  # noqa: E402
import _config as K  # noqa: E402
import rutirc as R  # noqa: E402

FAILS = []
T = pd.Timestamp


def check(name, cond, detail=""):
    print(("ok   " if cond else "FAIL ") + name + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        FAILS.append(name)


# =============================================================================
# PART 1. ALGORITHMS
# =============================================================================
def codes_df(rows):
    d = pd.DataFrame(rows, columns=["person_id", "date", "kind"])
    d["date"] = pd.to_datetime(d["date"])
    for k in ("upper", "lower", "unspec", "hematuria"):
        d[k] = d["kind"].eq(k)
    return d.drop(columns="kind")


def dates(rows):
    d = pd.DataFrame(rows, columns=["person_id", "date"]); d["date"] = pd.to_datetime(d["date"]); return d


def test_code_rx():
    m = R.match_codes(pd.Series(["N30.00", "N30.01", "N30.1", "N10", "N10.1", "C44.9", "C43.9", "Z85.828", "Z85.3"]),
                      ["N30.0*", "N10", "C4*", "Z85*"], ["C44*", "Z85.828"])
    check("code_rx: prefix, exact and exclusions", m.tolist() == [True, True, False, True, False, False, True, False, True], str(m.tolist()))


def test_episodes():
    codes = codes_df([(1, "2015-01-01", "lower"), (1, "2015-01-20", "upper"), (1, "2015-01-31", "lower"),
                      (1, "2015-02-01", "lower"), (1, "2015-06-01", "unspec"), (2, "2015-01-01", "lower")])
    abx = dates([(1, "2015-01-03"), (1, "2015-02-04")])
    cul = dates([(2, "2015-01-08")])
    E = R.build_episodes(codes, abx, cul, gap_days=31, abx_days=3, culture_days=7)
    e1 = E.loc[E["person_id"].eq(1)].sort_values("ep_start")
    check("episodes: fixed 31-day window (day 30 joins, day 31 opens a new episode)",
          e1["ep_start"].tolist() == [T("2015-01-01"), T("2015-02-01"), T("2015-06-01")], str(e1["ep_start"].tolist()))
    check("episodes: upper dominates within an episode", e1["ep_type"].tolist() == ["upper", "lower", "lower"])
    check("episodes: confirmation by antibiotic +/-3 d, unconfirmed episode not qualifying",
          e1["qualifying"].tolist() == [True, True, False], str(e1["qualifying"].tolist()))
    e2 = E.loc[E["person_id"].eq(2)]
    check("episodes: culture 7 days later confirms", bool(e2["qualifying"].iloc[0]) and bool(e2["confirmed_culture"].iloc[0]))


def test_recurrence():
    E = pd.DataFrame({"person_id": [1, 1, 2, 2, 2, 3, 3, 4, 5, 5],
                      "ep_start": pd.to_datetime(["2016-01-01", "2016-06-29", "2016-01-01", "2016-05-01", "2016-12-30",
                                                  "2016-01-01", "2016-08-01", "2016-01-01", "2014-06-01", "2015-03-01"]),
                      "ep_type": ["lower", "lower", "lower", "upper", "lower", "lower", "lower", "lower", "lower", "lower"],
                      "qualifying": [True] * 10})
    idx = R.index_candidates(E)
    check("index: an episode 273 d after another is not a clean index", not ((idx["person_id"] == 5) & (idx["E1"] == T("2015-03-01"))).any())
    check("index: the earlier one is", ((idx["person_id"] == 5) & (idx["E1"] == T("2014-06-01"))).any())
    r = R.recurrence_in_window(E, idx.loc[idx["person_id"].isin([1, 2, 3, 4])].reset_index(drop=True)).set_index("person_id")
    check("recurrence: 2 episodes 180 d apart = recurrent (O3 inclusive)", r.loc[1, "q1_group"] == "recurrent", r.loc[1, "q1_group"])
    check("recurrence: 3 episodes in 365 d = recurrent, any-upper", r.loc[2, "q1_group"] == "recurrent" and r.loc[2, "uti_type_window"] == "any_upper")
    check("recurrence: 2 episodes 213 d apart = intermediate", r.loc[3, "q1_group"] == "intermediate")
    check("recurrence: 1 episode = single", r.loc[4, "q1_group"] == "single")


def test_courses():
    rec = pd.DataFrame({"person_id": [1, 1, 1, 1, 2, 2], "abx_class": ["A", "A", "A", "C", "C", "C"],
                        "start_date": pd.to_datetime(["2016-01-01", "2016-01-10", "2016-02-01", "2016-01-05", "2016-03-01", "2016-03-02"]),
                        "rec_end": pd.to_datetime(["2016-01-06", "2016-01-15", "2016-02-06", "2016-01-12", "2016-03-02", "2016-03-03"]),
                        "visit_occurrence_id": [10, 11, 12, 13, 20, 20]})
    vis = pd.DataFrame({"visit_occurrence_id": [10, 11, 12, 13, 20], "visit_type": ["Outpatient Visit"] * 4 + ["Inpatient Visit"],
                        "visit_start_date": pd.to_datetime(["2016-01-01", "2016-01-10", "2016-02-01", "2016-01-05", "2016-03-01"]),
                        "visit_end_date": pd.to_datetime(["2016-01-01", "2016-01-10", "2016-02-01", "2016-01-05", "2016-03-05"])})
    codes = dates([(1, "2016-01-03")])
    Cc = R.build_courses(rec, vis, codes, merge_days=7, ind_days=3)
    a = Cc.loc[(Cc["person_id"] == 1) & (Cc["abx_class"] == "A")].sort_values("c_start")
    check("courses: start 4 d after previous end merges; 17 d after does not", len(a) == 2 and a["n_records"].tolist() == [2, 1], str(a))
    check("courses: urinary indication when a UTI code is within +/-3 d", a["urinary"].tolist() == [True, False])
    inp = Cc.loc[Cc["person_id"] == 2]
    check("courses: all records of a class in one inpatient stay = one course", len(inp) == 1 and inp["setting"].iloc[0] == "inpatient")


def test_prophylaxis():
    starts = pd.date_range("2017-01-01", periods=5, freq="30D")
    rec = pd.DataFrame({"person_id": 1, "prophylaxis_agent": "nitrofurantoin", "start_date": starts,
                        "rec_end": starts + pd.Timedelta(days=30)})
    P = R.prophylaxis_periods(rec, dates([(9, "2017-01-01")]), None)
    check("prophylaxis: 150 d continuous coverage detected", (P["rule"] == "coverage").any(), str(P))
    cs = pd.DataFrame({"person_id": 1, "abx_class": "A", "c_start": starts, "c_end": starts + pd.Timedelta(days=30),
                       "setting": "outpatient", "urinary": False, "n_records": 1})
    B = R.apply_prophylaxis_blocks(cs, P.loc[P["rule"].eq("coverage")])
    check("prophylaxis: courses replaced by one per 90-day block", len(B) == 2 and B["prophylaxis_block"].all(), str(len(B)))


def test_rcc():
    dx = pd.DataFrame([(1, "2018-01-01", "RCC"), (1, "2018-03-01", "RCC"),
                       (2, "2018-01-01", "RCC"),
                       (3, "2018-01-01", "RCC"), (3, "2018-01-10", "RENAL_PELVIS_CA"), (3, "2018-04-01", "RCC"),
                       (4, "2018-05-01", "RCC"), (4, "2018-04-20", "KIDNEY_DONOR"),
                       (5, "2018-05-01", "RCC"), (5, "2018-03-15", "RENAL_MASS"),
                       (6, "2018-01-01", "RCC"), (6, "2018-01-05", "RCC")], columns=["person_id", "date", "set_key"])
    dx["date"] = pd.to_datetime(dx["date"])
    px = pd.DataFrame([(4, "2018-05-10", "NEPHRECTOMY"), (5, "2018-04-25", "PARTIAL_NEPHRECTOMY")], columns=["person_id", "date", "set_key"])
    px["date"] = pd.to_datetime(px["date"])
    O = R.rcc_outcome(dx, px, dates([(9, "2018-01-01")])).set_index("person_id")
    check("rcc: second code 59 d later = confirmed", O.loc[1, "level"] == "confirmed" and O.loc[1, "confirm_route"] == "two_codes")
    check("rcc: one code only = possible", O.loc[2, "level"] == "possible")
    check("rcc: RCC + renal pelvis within 90 d, no procedure = ambiguous, not confirmed",
          O.loc[3, "ambiguous"] and O.loc[3, "level"] == "possible")
    check("rcc: nephrectomy next to a donor code does not confirm", O.loc[4, "level"] == "possible")
    check("rcc: partial nephrectomy confirms and dates the diagnosis",
          O.loc[5, "level"] == "confirmed" and O.loc[5, "dx_date"] == T("2018-04-25"), str(O.loc[5].to_dict()))
    check("rcc: detection date = renal mass code in the 180 d before", O.loc[5, "detection_date"] == T("2018-03-15"))
    check("rcc: two codes 4 d apart are not confirmation (rule-out)", O.loc[6, "level"] == "possible")


def test_egfr():
    e = R.egfr_ckdepi2021([1.0, 1.0], [60, 60], [False, True])
    check("eGFR CKD-EPI 2021: male 60 y, Scr 1.0 ~ 86; female ~ 64", abs(e[0] - 86) < 1.5 and abs(e[1] - 64) < 1.5, str(e))
    check("creatinine: umol/L converted, unknown unit dropped",
          np.allclose(R.creatinine_mgdl([88.4, 1.1, 1.0], ["micromole per liter", "milligram per deciliter", "x"])[:2], [1.0, 1.1])
          and np.isnan(R.creatinine_mgdl([1.0], ["x"])[0]))


def test_history():
    idx = pd.DataFrame({"person_id": [1, 2], "T0": pd.to_datetime(["2016-01-01", "2016-01-01"])})
    vd = dates([(1, "2013-06-01"), (1, "2014-06-01"), (1, "2015-06-01"), (2, "2012-06-01"), (2, "2015-06-01")])
    sp = pd.DataFrame({"person_id": [1, 2], "first_ehr": pd.to_datetime(["2012-01-01", "2012-01-01"])})
    ok = R.history_ok(idx, vd, sp)
    check("history: a visit in each of the 3 prior years required", ok.tolist() == [True, False], str(ok.tolist()))


# =============================================================================
# PART 2. SYNTHETIC CDR + END-TO-END RUN
# =============================================================================
def concrete(code):
    return code[:-1] + ("0" if not code[:-1].endswith(".") else "0") if code.endswith("*") else code


def build_cdr(seed=7, n=600):
    import duckdb
    rng = np.random.default_rng(seed)
    cid = iter(range(10_000_000, 99_000_000))
    concept, anc, rel = [], [], []

    def add(code, name, vocab, domain, cls="x", std=None):
        i = next(cid); concept.append((i, name, domain, vocab, cls, std, None, code)); return i
    code_id = {}
    for key, d in K.DX.items():
        word = d["rx"].split("|")[0]
        for vocab, k in (("ICD10CM", "icd10"), ("ICD9CM", "icd9")):
            for c in d.get(k, []):
                cc = concrete(c)
                if (vocab, cc) not in code_id:
                    i = add(cc, f"synthetic {word} {cc}", vocab, "Condition")
                    s = add("S" + cc, f"SNOMED {word} {cc}", "SNOMED", "Condition", std="S")
                    rel.append((i, s, "Maps to")); code_id[(vocab, cc)] = i
    for key, d in K.PX.items():
        word = d["rx"].split("|")[0]
        for vocab, k in (("CPT4", "cpt"), ("ICD10PCS", "pcs"), ("ICD9Proc", "icd9p")):
            for c in d.get(k, []):
                if (vocab, c) not in code_id:
                    code_id[(vocab, c)] = add(c, f"synthetic {word} {c}", vocab, "Procedure", std="S" if vocab == "CPT4" else None)
    for key, d in K.LABS.items():
        for c in d["loinc"]:
            code_id[("LOINC", c)] = add(c, {"CREATININE": "Creatinine [Mass/volume] in Serum", "HBA1C": "Hemoglobin A1c",
                                            "BMI": "Body mass index", "URINE_CULTURE_LAB": "Bacteria identified in Urine by Culture"}[key],
                                        "LOINC", "Measurement", std="S")
    # drugs
    ing = {}
    names = sorted({x for v in K.ABX_CLASSES.values() for x in v} | set(K.PROPHYLAXIS_ONLY) | set(K.NONABSORBED)
                   | set(K.NEG_CONTROL_RX) | set(K.RCC_DRUGS) | {"metformin", "lisinopril", "prednisone", "tacrolimus"})
    drug_id = {}
    for nm in names:
        i = add(nm, nm, "RxNorm", "Drug", "Ingredient", "S"); ing[nm] = i; anc.append((i, i, 0, 0))
        for form in ("100 MG Oral Tablet", "Ophthalmic Solution"):
            j = add(f"{nm} {form}", f"{nm} {form}", "RxNorm", "Drug", "Clinical Drug", "S"); anc.append((i, j, 1, 1))
            drug_id[(nm, form)] = j
    tmpsmx = add("tmpsmx", "sulfamethoxazole 800 MG / trimethoprim 160 MG Oral Tablet", "RxNorm", "Drug", "Clinical Drug", "S")
    anc += [(ing["sulfamethoxazole"], tmpsmx, 1, 1), (ing["trimethoprim"], tmpsmx, 1, 1)]
    for atc, nm in (("A10", "metformin"), ("C09", "lisinopril"), ("H02AB", "prednisone"), ("L04", "tacrolimus")):
        a = add(atc, {"A10": "DRUGS USED IN DIABETES", "C09": "AGENTS ACTING ON THE RENIN-ANGIOTENSIN SYSTEM",
                      "H02AB": "Glucocorticoids", "L04": "IMMUNOSUPPRESSANTS"}[atc], "ATC", "Drug", "ATC", "C")
        anc += [(a, ing[nm], 1, 1), (a, drug_id[(nm, "100 MG Oral Tablet")], 2, 2)]
    vis_c = {t: add(t, t, "Visit", "Visit", std="S") for t in ("Outpatient Visit", "Inpatient Visit", "Emergency Room Visit")}
    dtype = {t: add(t, t, "Type Concept", "Type Concept") for t in ("Prescription written", "Inpatient administration")}
    gend = {"MALE": add("M", "MALE", "Gender", "Gender"), "FEMALE": add("F", "FEMALE", "Gender", "Gender")}
    unit = {"mg/dL": add("mg/dL", "milligram per deciliter", "UCUM", "Unit")}
    ppi = {}
    for q, var in K.SURVEY_ITEMS.items():
        ppi[q] = add(q, {"race_eth": "Race: What Race Ethnicity", "sex_at_birth": "Biological Sex At Birth: Sex At Birth",
                         "education": "Education Level: Highest Grade", "income": "Income: Annual Income",
                         "smoke_100": "Smoking: 100 Cigs Lifetime", "zip3": "Street Address: PII ZIP"}[var], "PPI", "Observation", "Question", "S")
    for cw in K.CROSSWALKS.values():
        for a in cw:
            if a not in ppi:
                ppi[a] = add(a, a.split("_")[-1], "PPI", "Observation", "Answer", "S")
    ppi["Hispanic_PuertoRican"] = add("Hispanic_PuertoRican", "Puerto Rican", "PPI", "Observation", "Answer", "S")
    ppi["Smoke_Yes"] = add("Smoke_Yes", "Yes", "PPI", "Observation", "Answer", "S")
    consent = add("ConsentPII", "Consent PII", "PPI", "Observation", "Module", "S")
    anc.append((consent, consent, 0, 0))

    # people
    P, V, CO, DE, PO, ME, OB, DEATH = [], [], [], [], [], [], [], []
    vid = iter(range(1, 10 ** 8)); oid = iter(range(1, 10 ** 8))

    def dx(p, d, vocab, code):
        CO.append((next(oid), p, d, code_id[(vocab, code)], next(vid) if False else None))

    def rx(p, d, nm, form="100 MG Oral Tablet", days=7, typ="Prescription written", cid_=None):
        DE.append((next(oid), p, cid_ or drug_id[(nm, form)], d, d + pd.Timedelta(days=days), days, None, dtype[typ]))
    for p in range(1, n + 1):
        yob = int(rng.integers(1935, 1970)); male = rng.random() < 0.4
        P.append((p, yob, int(rng.integers(1, 13)), 15, None, gend["MALE" if male else "FEMALE"]))
        for d in pd.date_range("2006-01-15", "2024-06-15", freq="60D"):
            V.append((next(vid), p, d, d, vis_c["Outpatient Visit"]))
        OB.append((next(oid), p, consent, T("2018-05-01"), None, None))
        OB.append((next(oid), p, ppi["Race_WhatRaceEthnicity"], T("2018-05-01"),
                   ppi[rng.choice(["WhatRaceEthnicity_White", "WhatRaceEthnicity_Hispanic", "WhatRaceEthnicity_Black", "PMI_Skip"])], None))
        OB.append((next(oid), p, ppi["BiologicalSexAtBirth_SexAtBirth"], T("2018-05-01"),
                   ppi["SexAtBirth_Male" if male else "SexAtBirth_Female"], None))
        OB.append((next(oid), p, ppi["StreetAddress_PIIZIP"], T("2018-05-01"), None, rng.choice(["00612", "10001", "00901"])))
        OB.append((next(oid), p, ppi["Smoking_100CigsLifetime"], T("2018-05-01"), ppi["Smoke_Yes"] if rng.random() < .3 else None, None))
        # UTI history: index episode 2011-2016, 0-3 further episodes in the year after
        e1 = T("2011-01-01") + pd.Timedelta(days=int(rng.integers(0, 2000)))
        k = int(rng.choice([0, 1, 2, 3], p=[.5, .2, .15, .15]))
        eps = [e1] + [e1 + pd.Timedelta(days=int(x)) for x in sorted(rng.choice(np.arange(40, 360, 40), k, replace=False))]
        for j, d in enumerate(eps):
            icd10 = d >= T(K.ICD10_START)
            upper = rng.random() < .15
            code = ("N10" if upper else rng.choice(["N30.00", "N39.0"])) if icd10 else ("590.10" if upper else "599.0")
            dx(p, d, "ICD10CM" if icd10 else "ICD9CM", code)
            if rng.random() < .85:
                rx(p, d + pd.Timedelta(days=int(rng.integers(0, 3))), str(rng.choice(["nitrofurantoin", "ciprofloxacin", "cephalexin"])), days=5)
            elif rng.random() < .5:
                rx(p, d, "", cid_=tmpsmx, days=3)
        if rng.random() < .1:
            rx(p, e1 + pd.Timedelta(days=500), "nitrofurantoin", form="Ophthalmic Solution")
        if rng.random() < .3:
            rx(p, e1 - pd.Timedelta(days=400), "metformin", days=90); rx(p, e1 - pd.Timedelta(days=300), "lisinopril", days=90)
            rx(p, e1 - pd.Timedelta(days=200), "lisinopril", days=90)
        ME.append((next(oid), p, code_id[("LOINC", "2160-0")], e1 - pd.Timedelta(days=100), float(rng.uniform(.6, 1.6)), unit["mg/dL"], None))
        ME.append((next(oid), p, code_id[("LOINC", "39156-5")], e1 - pd.Timedelta(days=200), float(rng.uniform(19, 42)), None, None))
        if rng.random() < .08:                                     # incident RCC, confirmed by a second code or surgery
            d = e1 + pd.Timedelta(days=int(rng.integers(900, 3000)))
            if d < T("2024-12-01"):
                dx(p, d, "ICD10CM" if d >= T(K.ICD10_START) else "ICD9CM", "C640" if d >= T(K.ICD10_START) else "189.0")
                if rng.random() < .7:
                    dx(p, d + pd.Timedelta(days=60), "ICD10CM", "C640")
                else:
                    PO.append((next(oid), p, d + pd.Timedelta(days=20), code_id[("CPT4", "50240")], code_id[("CPT4", "50240")]))
        if rng.random() < .03:
            dx(p, e1 - pd.Timedelta(days=50), "ICD10CM" if e1 >= T(K.ICD10_START) else "ICD9CM",
               "C50*".replace("*", "0") if e1 >= T(K.ICD10_START) else "174*".replace("*", "0"))
        if rng.random() < .05:
            DEATH.append((p, e1 + pd.Timedelta(days=int(rng.integers(800, 4000))), True))
        if rng.random() < .1:
            PO.append((next(oid), p, e1 - pd.Timedelta(days=300), code_id[("CPT4", "74176")], code_id[("CPT4", "74176")]))
    # CO visit ids, ext tables
    db = duckdb.connect(); db.execute('CREATE SCHEMA "cdr"')
    frames = {
        "concept": pd.DataFrame(concept, columns=["concept_id", "concept_name", "domain_id", "vocabulary_id", "concept_class_id",
                                                  "standard_concept", "invalid_reason", "concept_code"]),
        "concept_ancestor": pd.DataFrame(anc, columns=["ancestor_concept_id", "descendant_concept_id", "min_levels_of_separation", "max_levels_of_separation"]),
        "concept_relationship": pd.DataFrame(rel, columns=["concept_id_1", "concept_id_2", "relationship_id"]),
        "person": pd.DataFrame(P, columns=["person_id", "year_of_birth", "month_of_birth", "day_of_birth", "birth_datetime", "gender_concept_id"]),
        "visit_occurrence": pd.DataFrame(V, columns=["visit_occurrence_id", "person_id", "visit_start_date", "visit_end_date", "visit_concept_id"]),
        "condition_occurrence": pd.DataFrame(CO, columns=["condition_occurrence_id", "person_id", "condition_start_date", "condition_source_concept_id", "visit_occurrence_id"]),
        "drug_exposure": pd.DataFrame(DE, columns=["drug_exposure_id", "person_id", "drug_concept_id", "drug_exposure_start_date",
                                                   "drug_exposure_end_date", "days_supply", "visit_occurrence_id", "drug_type_concept_id"]),
        "procedure_occurrence": pd.DataFrame(PO, columns=["procedure_occurrence_id", "person_id", "procedure_date", "procedure_source_concept_id", "procedure_concept_id"]),
        "measurement": pd.DataFrame(ME, columns=["measurement_id", "person_id", "measurement_concept_id", "measurement_date", "value_as_number", "unit_concept_id", "value_as_concept_id"]),
        "observation": pd.DataFrame(OB, columns=["observation_id", "person_id", "observation_source_concept_id", "observation_date", "value_source_concept_id", "value_as_string"]),
        "aou_death": pd.DataFrame(DEATH, columns=["person_id", "death_date", "primary_death_record"]),
    }
    f = frames
    f["drug_exposure"]["route_concept_id"] = 0; f["drug_exposure"]["quantity"] = None
    f["procedure_occurrence"]["visit_occurrence_id"] = None
    f["measurement"]["measurement_source_concept_id"] = f["measurement"]["measurement_concept_id"]; f["measurement"]["visit_occurrence_id"] = None
    f["observation"]["visit_occurrence_id"] = None
    f["observation"]["observation_concept_id"] = f["observation"]["observation_source_concept_id"]
    f["condition_occurrence"]["visit_occurrence_id"] = f["condition_occurrence"]["visit_occurrence_id"].astype("Int64")
    f["drug_exposure"]["visit_occurrence_id"] = f["drug_exposure"]["visit_occurrence_id"].astype("Int64")
    f["aou_death"]["src_id"] = "EHR site 1"
    f["observation_period"] = f["person"][["person_id"]].assign(observation_period_id=1)
    for t, k in (("visit_occurrence", "visit_occurrence_id"), ("condition_occurrence", "condition_occurrence_id"),
                 ("drug_exposure", "drug_exposure_id"), ("procedure_occurrence", "procedure_occurrence_id"),
                 ("measurement", "measurement_id"), ("observation", "observation_id")):
        f[f"{t}_ext"] = pd.DataFrame({k: f[t][k], "src_id": np.where(f[t]["person_id"] % 2 == 0, "EHR site 101", "EHR site 202")})
    for c in ("visit_start_date", "visit_end_date"):
        f["visit_occurrence"][c] = pd.to_datetime(f["visit_occurrence"][c]).dt.date
    for t, c in (("condition_occurrence", "condition_start_date"), ("drug_exposure", "drug_exposure_start_date"),
                 ("drug_exposure", "drug_exposure_end_date"), ("procedure_occurrence", "procedure_date"),
                 ("measurement", "measurement_date"), ("observation", "observation_date"), ("aou_death", "death_date")):
        f[t][c] = pd.to_datetime(f[t][c]).dt.date
    for nme, d in f.items():
        db.register("_t", d); db.execute(f'CREATE TABLE "cdr"."{nme}" AS SELECT * FROM _t'); db.unregister("_t")
    return db


def make_q(db):
    import sqlglot

    def lit(v):
        if isinstance(v, str):
            return "'" + v.replace("\\", "\\\\").replace("'", "\\'") + "'"
        return repr(v)

    def q_override(sql, prms, label):
        s = sql
        for p in sorted(prms, key=lambda p: -len(p.name)):
            s = s.replace("@" + p.name, ("[" + ", ".join(lit(v) for v in p.value) + "]") if p.array else lit(p.value))
        s = sqlglot.transpile(s, read="bigquery", write="duckdb")[0]
        return db.execute(s).df()
    return q_override


def scan_tables(run):
    """Every saved CSV: no count column holds 1..20 (all count-like columns per aou_core.COUNT_RX)."""
    bad = []
    for f in sorted((run / "tables").glob("*.csv")):
        if f.name == "_index.csv":
            continue
        d = pd.read_csv(f)
        for c in d.columns:
            if re.search(C.COUNT_RX, c):
                v = pd.to_numeric(d[c], errors="coerce")
                if ((v >= 1) & (v <= 20)).any():
                    bad.append(f"{f.name}:{c}")
    return bad


def test_end_to_end():
    db = build_cdr()
    tmp = Path(tempfile.mkdtemp())
    C.CFG.update(run_dir=str(tmp / "run"), q_override=make_q(db), cdr="cdr", version="TEST")
    RUN_CONTROL = {"new_run": True}
    code_hash = hashlib.sha256(b"test").hexdigest()[:12]
    for cell in ["F00_setup", "F01_inventory", "F02_extract", "F03_phenotypes", "F04_cohorts"]:
        C.STATE["failed"] = []
        exec(compile((HERE / f"{cell}.py").read_text(), cell, "exec"), {"RUN_CONTROL": RUN_CONTROL, "CODE_HASH": code_hash,
                                                                        "__name__": "__main__"})
        check(f"end-to-end: {cell} ran with no failed part", not C.STATE["failed"], "; ".join(C.STATE["failed"]))
    run = Path(C.CFG["run_dir"])
    F1 = C.load_state("frame_Q1")
    check("end-to-end: Q1 cohort has recurrent and single persons", set(F1["q1_group"]) >= {"recurrent", "single"})
    check("end-to-end: Q1 T0 = E1 + 365 d", ((F1["T0"] - F1["E1"]).dt.days == 365).all())
    check("end-to-end: outcome counting starts 2 y after T0; no event before it",
          (F1["S"] > F1["T0"]).all() and not (F1["dx_date"] < F1["S"]).any())
    check("end-to-end: follow-up ends by the CDR cutoff", (F1["end"] <= T(K.CDR_CUTOFF)).all())
    check("end-to-end: age >= 40 at T0", (F1["age_t0"] >= 40).all())
    F2 = C.load_state("frame_Q2")
    check("end-to-end: Q2 classes in {A, B, C, D}", set(F2["q2_class"]) <= {"A", "B", "C", "D"} and len(F2) > 0)
    G = C.load_state("drug_glossary", internal=False)
    check("end-to-end: ophthalmic antibiotic not counted", not G.loc[G["form"].eq("non_systemic"), "abx_counted"].any())
    check("end-to-end: TMP-SMX product is class B once", G.loc[G["concept_name"].str.contains("sulfamethoxazole 800"), "abx_classes"].iloc[0] == "B")
    bad = scan_tables(run)
    check("disclosure: no saved table shows a count 1-20", not bad, ", ".join(bad[:10]))
    txt = (run / "tables" / "ALL_TABLES.txt").read_text()
    check("disclosure: ALL_TABLES.txt has no person_id column", "person_id" not in txt)
    check("blind: no event count by exposure group in F04 tables",
          not any("event" in c for f in (run / "tables").glob("F04_table1.csv") for c in pd.read_csv(f).columns))


if __name__ == "__main__":
    for t in [test_code_rx, test_episodes, test_recurrence, test_courses, test_prophylaxis, test_rcc, test_egfr, test_history,
              test_end_to_end]:
        print(f"\n## {t.__name__}")
        t()
    print("\n" + ("all tests passed" if not FAILS else f"{len(FAILS)} FAILED: " + "; ".join(FAILS)))
    sys.exit(1 if FAILS else 0)
