"""rutirc.py - RUTI-RCC study functions on top of aou_core (code-list resolution, SQL, phenotype algorithms).

Everything definitional reads _config.py (imported as K). Pure pandas functions below the SQL section are the
protocol's algorithms; each cites the section it implements and is exercised by test_synthetic.py.
No participant-level data is ever written outside the run folder (aou_core.save_state names it _INTERNAL).
"""
import re

import numpy as np
import pandas as pd

import aou_core as C
import _config as K

DAY = pd.Timedelta(days=1)


# =============================================================================
# 1. CODE LISTS -> REGEX, AND RESOLUTION INSIDE THE CDR (F00)
# =============================================================================
def code_rx(codes, exclude=()):
    """["N30.0*", "N10"] -> anchored regex. '*' = the code and everything under it (prefix)."""
    def one(c):
        c = c.strip()
        return "^" + re.escape(c[:-1]) if c.endswith("*") else "^" + re.escape(c) + "$"
    inc = "|".join(one(c) for c in codes) if codes else "a^"
    exc = "|".join(one(c) for c in exclude) if exclude else None
    return inc, exc


def match_codes(series, codes, exclude=()):
    inc, exc = code_rx(codes, exclude)
    s = series.astype(str)
    m = s.str.contains(inc, regex=True)
    if exc:
        m &= ~s.str.contains(exc, regex=True)
    return m


CODE_VOCABS = ["ICD9CM", "ICD10CM", "CPT4", "ICD10PCS", "ICD9Proc", "LOINC", "SNOMED", "ATC"]
CONCEPT_BY_CODE_SQL = """SELECT concept_id, concept_code, concept_name, vocabulary_id, domain_id, concept_class_id,
       standard_concept, invalid_reason
FROM `{CDR}.concept`
WHERE vocabulary_id IN UNNEST(@vocabs) AND REGEXP_CONTAINS(concept_code, @rx)"""
MAPS_TO_SQL = """SELECT r.concept_id_1 AS source_concept_id, c.concept_id AS standard_concept_id,
       c.concept_name AS standard_name, c.domain_id AS standard_domain, c.vocabulary_id AS standard_vocabulary
FROM `{CDR}.concept_relationship` r JOIN `{CDR}.concept` c ON c.concept_id = r.concept_id_2
WHERE r.relationship_id = 'Maps to' AND r.concept_id_1 IN UNNEST(@ids)"""
INGREDIENT_SQL = """SELECT concept_id, concept_name, concept_class_id, standard_concept, vocabulary_id
FROM `{CDR}.concept`
WHERE vocabulary_id = 'RxNorm' AND concept_class_id IN ('Ingredient', 'Precise Ingredient')
  AND LOWER(concept_name) IN UNNEST(@names)"""
DRUG_DESC_SQL = """SELECT ca.ancestor_concept_id AS root_id, c.concept_id, c.concept_name, c.concept_class_id,
       c.vocabulary_id, c.standard_concept
FROM `{CDR}.concept_ancestor` ca JOIN `{CDR}.concept` c ON c.concept_id = ca.descendant_concept_id
WHERE ca.ancestor_concept_id IN UNNEST(@ids) AND c.domain_id = 'Drug'"""


def all_code_regex():
    """One regex over every code the study names, so F00 resolves all code sets in one concept-table scan."""
    pats = []
    for d in K.DX.values():
        for k in ("icd10", "icd9"):
            pats.append(code_rx(d.get(k, []))[0])
    for d in K.PX.values():
        for k in ("cpt", "pcs", "icd9p"):
            if d.get(k):
                pats.append(code_rx(d[k])[0])
    for d in K.LABS.values():
        pats.append(code_rx(d["loinc"] + d.get("snomed", []))[0])
    for d in K.ATC_GROUPS.values():
        pats.append(code_rx(d["atc"])[0])
    return "|".join(f"(?:{p})" for p in pats if p != "a^")


def assign_code_sets(con):
    """concept rows (one scan) -> long table set_key, kind, concept rows. ICD-9 and ICD-10 kept apart (V/E clash)."""
    out = []
    for key, d in K.DX.items():
        for vocab, k in (("ICD10CM", "icd10"), ("ICD9CM", "icd9")):
            sub = con.loc[con["vocabulary_id"].eq(vocab)]
            m = match_codes(sub["concept_code"], d.get(k, []), d.get(f"{k}_exclude", ()))
            out.append(sub.loc[m].assign(set_key=key, kind="DX", rx_expected=d["rx"]))
    for key, d in K.PX.items():
        for vocab, k in (("CPT4", "cpt"), ("ICD10PCS", "pcs"), ("ICD9Proc", "icd9p")):
            sub = con.loc[con["vocabulary_id"].eq(vocab)]
            m = match_codes(sub["concept_code"], d.get(k, []))
            out.append(sub.loc[m].assign(set_key=key, kind="PX", rx_expected=d["rx"]))
    for key, d in K.LABS.items():
        for vocab, k in (("LOINC", "loinc"), ("SNOMED", "snomed")):
            sub = con.loc[con["vocabulary_id"].eq(vocab)]
            out.append(sub.loc[match_codes(sub["concept_code"], d.get(k, []))].assign(set_key=key, kind="LAB",
                                                                                         rx_expected=d["rx"]))
    for key, d in K.ATC_GROUPS.items():
        sub = con.loc[con["vocabulary_id"].eq("ATC")]
        out.append(sub.loc[match_codes(sub["concept_code"], d["atc"])].assign(set_key=key, kind="ATC", rx_expected=d["rx"]))
    S = pd.concat(out, ignore_index=True)
    S["name_ok"] = [bool(re.search(r, str(n), re.I)) for r, n in zip(S["rx_expected"], S["concept_name"])]
    S["era"] = np.where(S["vocabulary_id"].isin(["ICD9CM", "ICD9Proc"]), "ICD-9",
                        np.where(S["vocabulary_id"].isin(["ICD10CM", "ICD10PCS"]), "ICD-10", "n/a"))
    return S


def code_set_coverage(S):
    """Every code the protocol names, and whether the CDR has it (codes absent from the CDR are listed, not dropped
    silently). Wildcards count as found when at least one code under them exists."""
    rows = []
    for key, d in list(K.DX.items()) + list(K.PX.items()):
        for k, vocab in (("icd10", "ICD10CM"), ("icd9", "ICD9CM"), ("cpt", "CPT4"), ("pcs", "ICD10PCS"), ("icd9p", "ICD9Proc")):
            for code in d.get(k, []):
                sub = S.loc[S["set_key"].eq(key) & S["vocabulary_id"].eq(vocab)]
                hit = match_codes(sub["concept_code"], [code]).sum() if len(sub) else 0
                rows.append({"set_key": key, "vocabulary": vocab, "code": code, "codes_found": int(hit)})
    return pd.DataFrame(rows)


def classify_drug_form(name):
    n = str(name).lower()
    if re.search(K.NON_SYSTEMIC_FORM_RX, n):
        return "non_systemic"
    if re.search(r"inject|intraven|intramusc|\biv\b|infusion|for solution|vial", n):
        return "parenteral"
    if re.search(K.ORAL_FORM_RX, n):
        return "oral"
    return "unknown_form"


def drug_glossary(ing, desc):
    """Ingredients (by name) + descendants -> one row per drug concept with: ingredient list, study groups (ABX class
    letters, NEG_CONTROL, RCC_DRUG, PROPHYLAXIS_ONLY, NONABSORBED, URINARY_ACTIVE flag), form and the counted flag."""
    ing = ing.copy()
    ing["name_l"] = ing["concept_name"].str.lower()
    id2name = dict(zip(ing["concept_id"].astype(int), ing["name_l"]))
    d = desc.copy()
    d["ingredient"] = d["root_id"].astype(int).map(id2name)
    g = d.groupby(["concept_id", "concept_name", "concept_class_id"], as_index=False)["ingredient"].agg(
        lambda s: sorted(set(s)))
    cls_of = {n.lower(): c for c, names in K.ABX_CLASSES.items() for n in names}
    ua = {n.lower() for n in K.URINARY_ACTIVE}

    def classes(ings):
        return sorted({cls_of[i] for i in ings if i in cls_of and i not in K.COMPANION_ONLY}) or \
            sorted({cls_of[i] for i in ings if i in cls_of})
    g["abx_classes"] = g["ingredient"].map(classes)
    g["form"] = g["concept_name"].map(classify_drug_form)
    g["is_abx"] = g["abx_classes"].map(len).gt(0)
    oral_vanco = g["ingredient"].map(lambda s: "vancomycin" in s) & g["concept_name"].str.contains(K.ORAL_FORM_RX, case=False, regex=True)
    nonabs = g["ingredient"].map(lambda s: bool(set(s) & set(K.NONABSORBED))) | oral_vanco
    g["abx_counted"] = g["is_abx"] & g["form"].ne("non_systemic") & ~nonabs
    g["urinary_active"] = g["abx_counted"] & g["ingredient"].map(lambda s: bool(set(s) & ua))
    g["neg_control"] = g["ingredient"].map(lambda s: bool(set(s) & set(K.NEG_CONTROL_RX))) & g["form"].ne("non_systemic")
    g["rcc_drug"] = g["ingredient"].map(lambda s: bool(set(s) & set(K.RCC_DRUGS)))
    g["prophylaxis_only"] = g["ingredient"].map(lambda s: bool(set(s) & set(K.PROPHYLAXIS_ONLY)))
    g["prophylaxis_agent"] = g["ingredient"].map(lambda s: next((a for a in K.PROPHYLAXIS_AGENTS if a in s), None))
    g["multi_class"] = g["abx_classes"].map(len).gt(1)
    g["ingredient"] = g["ingredient"].map(lambda s: "+".join(s))
    g["abx_classes"] = g["abx_classes"].map(lambda s: "+".join(s))
    return g


# =============================================================================
# 2. EXTRACTION SQL (F02). {CDR} is filled by aou_core.q; @ids are resolved concept ids from F00 state.
# =============================================================================
CAND_CTE = """cand AS (
  SELECT DISTINCT person_id FROM `{CDR}.condition_occurrence` WHERE condition_source_concept_id IN UNNEST(@cand_ids)
  UNION DISTINCT
  SELECT DISTINCT person_id FROM `{CDR}.observation` WHERE observation_source_concept_id IN UNNEST(@cand_ids)
)"""

DX_EVENTS_SQL = "WITH " + CAND_CTE + """
SELECT x.person_id, x.event_date, x.source_concept_id, x.visit_occurrence_id, x.tbl
FROM (
  SELECT person_id, condition_start_date AS event_date, condition_source_concept_id AS source_concept_id,
         visit_occurrence_id, 'condition' AS tbl
  FROM `{CDR}.condition_occurrence` WHERE condition_source_concept_id IN UNNEST(@ids)
  UNION ALL
  SELECT person_id, observation_date, observation_source_concept_id, visit_occurrence_id, 'observation'
  FROM `{CDR}.observation` WHERE observation_source_concept_id IN UNNEST(@ids)
  UNION ALL
  SELECT person_id, procedure_date, procedure_source_concept_id, visit_occurrence_id, 'procedure'
  FROM `{CDR}.procedure_occurrence` WHERE procedure_source_concept_id IN UNNEST(@ids)
  UNION ALL
  SELECT person_id, measurement_date, measurement_source_concept_id, visit_occurrence_id, 'measurement'
  FROM `{CDR}.measurement` WHERE measurement_source_concept_id IN UNNEST(@ids)
) x JOIN cand USING (person_id)"""

PX_EVENTS_SQL = "WITH " + CAND_CTE + """
SELECT x.person_id, x.event_date, x.concept_id, x.visit_occurrence_id, x.tbl
FROM (
  SELECT person_id, procedure_date AS event_date, procedure_source_concept_id AS concept_id, visit_occurrence_id,
         'procedure_src' AS tbl
  FROM `{CDR}.procedure_occurrence` WHERE procedure_source_concept_id IN UNNEST(@ids)
  UNION ALL
  SELECT person_id, procedure_date, procedure_concept_id, visit_occurrence_id, 'procedure_std'
  FROM `{CDR}.procedure_occurrence` WHERE procedure_concept_id IN UNNEST(@ids)
  UNION ALL
  SELECT person_id, measurement_date, measurement_source_concept_id, visit_occurrence_id, 'measurement'
  FROM `{CDR}.measurement` WHERE measurement_source_concept_id IN UNNEST(@ids)
  UNION ALL
  SELECT person_id, observation_date, observation_source_concept_id, visit_occurrence_id, 'observation'
  FROM `{CDR}.observation` WHERE observation_source_concept_id IN UNNEST(@ids)
) x JOIN cand USING (person_id)"""

DRUG_SQL = "WITH " + CAND_CTE + """
SELECT d.person_id, d.drug_concept_id, d.drug_exposure_start_date AS start_date, d.drug_exposure_end_date AS end_date,
       d.days_supply, d.quantity, d.visit_occurrence_id, d.route_concept_id, t.concept_name AS drug_type
FROM `{CDR}.drug_exposure` d JOIN cand USING (person_id)
LEFT JOIN `{CDR}.concept` t ON t.concept_id = d.drug_type_concept_id
WHERE d.drug_concept_id IN UNNEST(@ids)"""

VISITS_SQL = "WITH " + CAND_CTE + """
SELECT v.person_id, v.visit_occurrence_id, v.visit_start_date, v.visit_end_date, c.concept_name AS visit_type, e.src_id
FROM `{CDR}.visit_occurrence` v JOIN cand USING (person_id)
LEFT JOIN `{CDR}.concept` c ON c.concept_id = v.visit_concept_id
LEFT JOIN `{CDR}.visit_occurrence_ext` e ON e.visit_occurrence_id = v.visit_occurrence_id"""

# P§3.2: EHR activity = any EHR-sourced record in visit, condition, procedure, drug, measurement or observation.
_SPAN_PART = """  SELECT t.person_id, MIN(t.{d}) AS first_d, MAX(t.{d}) AS last_d
  FROM `{{CDR}}.{tbl}` t JOIN cand USING (person_id)
  JOIN `{{CDR}}.{tbl}_ext` e ON e.{tbl}_id = t.{tbl}_id
  WHERE REGEXP_CONTAINS(COALESCE(e.src_id, ''), @ehr_rx)
  GROUP BY t.person_id"""
EHR_SPAN_SQL = "WITH " + CAND_CTE + ", parts AS (\n" + "\n  UNION ALL\n".join(
    _SPAN_PART.format(tbl=t, d=d) for t, d in [
        ("visit_occurrence", "visit_start_date"), ("condition_occurrence", "condition_start_date"),
        ("procedure_occurrence", "procedure_date"), ("drug_exposure", "drug_exposure_start_date"),
        ("measurement", "measurement_date"), ("observation", "observation_date")]) + """
)
SELECT person_id, MIN(first_d) AS first_ehr, MAX(last_d) AS last_ehr FROM parts GROUP BY person_id"""

LABS_SQL = "WITH " + CAND_CTE + """
SELECT m.person_id, m.measurement_date, m.measurement_concept_id, m.value_as_number, u.concept_name AS unit,
       v.concept_name AS value_as_concept
FROM `{CDR}.measurement` m JOIN cand USING (person_id)
LEFT JOIN `{CDR}.concept` u ON u.concept_id = m.unit_concept_id
LEFT JOIN `{CDR}.concept` v ON v.concept_id = m.value_as_concept_id
WHERE m.measurement_concept_id IN UNNEST(@ids)"""

PERSON_SQL = "WITH " + CAND_CTE + """
SELECT p.person_id, p.year_of_birth, p.month_of_birth, p.day_of_birth, p.birth_datetime, g.concept_name AS gender
FROM `{CDR}.person` p JOIN cand USING (person_id)
LEFT JOIN `{CDR}.concept` g ON g.concept_id = p.gender_concept_id"""

SURVEY_SQL = "WITH " + CAND_CTE + """
SELECT o.person_id, o.observation_date, q.concept_code AS question_code, a.concept_code AS answer_code,
       a.concept_name AS answer_name, o.value_as_string
FROM `{CDR}.observation` o JOIN cand USING (person_id)
JOIN `{CDR}.concept` q ON q.concept_id = o.observation_source_concept_id
LEFT JOIN `{CDR}.concept` a ON a.concept_id = o.value_source_concept_id
WHERE q.concept_code IN UNNEST(@qcodes)
  AND COALESCE(o.value_source_concept_id, 0) != 2000000010"""

# F01 inventory (whole CDR, aggregate only)
SURVEY_INVENTORY_SQL = """SELECT q.concept_code AS question_code, q.concept_name AS question, a.concept_code AS answer_code,
       a.concept_name AS answer, COUNT(DISTINCT o.person_id) AS persons
FROM `{CDR}.observation` o
JOIN `{CDR}.concept` q ON q.concept_id = o.observation_source_concept_id
LEFT JOIN `{CDR}.concept` a ON a.concept_id = o.value_source_concept_id
WHERE q.vocabulary_id = 'PPI' AND REGEXP_CONTAINS(LOWER(q.concept_name), @rx)
GROUP BY 1, 2, 3, 4"""
DOB_SQL = """SELECT COUNT(*) AS persons, COUNTIF(day_of_birth IS NULL) AS n_day_missing,
       COUNTIF(day_of_birth = 15) AS n_day_15, COUNTIF(birth_datetime IS NULL) AS n_birth_datetime_missing,
       COUNT(DISTINCT day_of_birth) AS distinct_days,
       COUNT(DISTINCT EXTRACT(DAY FROM CAST(birth_datetime AS TIMESTAMP))) AS bdt_distinct_days,
       COUNT(DISTINCT EXTRACT(MONTH FROM CAST(birth_datetime AS TIMESTAMP))) AS bdt_distinct_months,
       AVG(IF(EXTRACT(DAY FROM CAST(birth_datetime AS TIMESTAMP)) = 1 AND EXTRACT(MONTH FROM CAST(birth_datetime AS TIMESTAMP)) = 1, 1, 0)) AS bdt_share_jan1,
       AVG(IF(EXTRACT(DAY FROM CAST(birth_datetime AS TIMESTAMP)) = 15, 1, 0)) AS bdt_share_day15
FROM `{CDR}.person`"""
SRC_PATTERN_SQL = """SELECT REGEXP_REPLACE(COALESCE(e.src_id, '(null)'), r'[0-9]+', '#') AS src_pattern,
       COUNT(DISTINCT e.src_id) AS distinct_src_ids, COUNT(DISTINCT v.person_id) AS persons
FROM `{CDR}.visit_occurrence` v LEFT JOIN `{CDR}.visit_occurrence_ext` e ON e.visit_occurrence_id = v.visit_occurrence_id
GROUP BY 1"""
DRUG_QUALITY_SQL = """SELECT t.concept_name AS drug_type, COUNT(*) AS n_records, COUNT(DISTINCT d.person_id) AS persons,
       AVG(IF(d.days_supply > 0, 1, 0)) AS share_days_supply, AVG(IF(d.drug_exposure_end_date IS NOT NULL, 1, 0)) AS share_end_date,
       AVG(IF(COALESCE(d.route_concept_id, 0) != 0, 1, 0)) AS share_route, AVG(IF(d.visit_occurrence_id IS NOT NULL, 1, 0)) AS share_visit
FROM `{CDR}.drug_exposure` d LEFT JOIN `{CDR}.concept` t ON t.concept_id = d.drug_type_concept_id
WHERE d.drug_concept_id IN UNNEST(@ids)
GROUP BY 1"""
DEATH_SOURCE_SQL = """SELECT REGEXP_REPLACE(COALESCE(CAST(src_id AS STRING), '(null)'), r'[0-9]+', '#') AS src_pattern,
       primary_death_record, COUNT(DISTINCT person_id) AS persons
FROM `{CDR}.aou_death` GROUP BY 1, 2"""
LAB_NAME_SEARCH_SQL = """SELECT c.concept_id, c.concept_code, c.concept_name, c.vocabulary_id,
       COUNT(DISTINCT m.person_id) AS persons, AVG(IF(m.value_as_concept_id IS NOT NULL AND m.value_as_concept_id != 0, 1, 0)) AS share_value_concept,
       AVG(IF(m.value_as_number IS NOT NULL, 1, 0)) AS share_value_number
FROM `{CDR}.measurement` m JOIN `{CDR}.concept` c ON c.concept_id = m.measurement_concept_id
WHERE REGEXP_CONTAINS(LOWER(c.concept_name), @rx)
GROUP BY 1, 2, 3, 4"""
SITE_YEAR_SQL = """WITH v AS (
  SELECT e.src_id, EXTRACT(YEAR FROM t.visit_start_date) AS yr, COUNT(*) AS n_visits
  FROM `{CDR}.visit_occurrence` t JOIN `{CDR}.visit_occurrence_ext` e ON e.visit_occurrence_id = t.visit_occurrence_id
  WHERE REGEXP_CONTAINS(COALESCE(e.src_id, ''), @ehr_rx) GROUP BY 1, 2),
d AS (
  SELECT e.src_id, EXTRACT(YEAR FROM t.drug_exposure_start_date) AS yr, COUNT(*) AS n_drugs
  FROM `{CDR}.drug_exposure` t JOIN `{CDR}.drug_exposure_ext` e ON e.drug_exposure_id = t.drug_exposure_id
  WHERE REGEXP_CONTAINS(COALESCE(e.src_id, ''), @ehr_rx) GROUP BY 1, 2)
SELECT v.src_id, v.yr, v.n_visits, COALESCE(d.n_drugs, 0) AS n_drugs
FROM v LEFT JOIN d ON d.src_id = v.src_id AND d.yr = v.yr"""
ZIP3_SQL = """SELECT SUBSTR(o.value_as_string, 1, 3) AS zip3, COUNT(DISTINCT o.person_id) AS persons
FROM `{CDR}.observation` o JOIN `{CDR}.concept` q ON q.concept_id = o.observation_source_concept_id
WHERE q.concept_code = @zip_code AND SUBSTR(o.value_as_string, 1, 3) IN UNNEST(@zips)
GROUP BY 1"""
OBS_PERIOD_SQL = """SELECT COUNT(DISTINCT person_id) AS persons, COUNT(*) AS n_periods FROM `{CDR}.observation_period`"""


# =============================================================================
# 3. SMALL HELPERS
# =============================================================================
def birth_dates(P):
    """P§3.1: full date of birth when released, else day 15 assumed (O6). Returns Series indexed by person_id."""
    bd = C.to_date(P["birth_datetime"]).values
    y = P["year_of_birth"].astype("Int64")
    m = P["month_of_birth"].fillna(7).astype("Int64") if "month_of_birth" in P else 7
    dd = P["day_of_birth"].fillna(K.IMPLEMENTATION["O6_dob_day_if_missing"]).astype("Int64")
    built = pd.to_datetime(pd.DataFrame({"year": y, "month": m, "day": dd}), errors="coerce")
    out = pd.Series(np.where(pd.isna(bd), built, bd), index=P["person_id"].values)
    return pd.to_datetime(out)


def egfr_ckdepi2021(scr_mgdl, age, female):
    """CKD-EPI 2021 creatinine equation (race-free). scr in mg/dL, age in years, female boolean."""
    scr = np.asarray(scr_mgdl, float); age = np.asarray(age, float); f = np.asarray(female, bool)
    kappa = np.where(f, 0.7, 0.9); alpha = np.where(f, -0.241, -0.302)
    r = scr / kappa
    return 142 * np.minimum(r, 1) ** alpha * np.maximum(r, 1) ** -1.200 * 0.9938 ** age * np.where(f, 1.012, 1.0)


def creatinine_mgdl(value, unit):
    """Units: mg/dL kept; umol/L divided by 88.4; anything else -> NaN (reported in F03)."""
    u = pd.Series(unit, dtype=object).astype(str).str.lower()
    v = pd.to_numeric(pd.Series(value), errors="coerce").values
    out = np.where(u.str.contains(r"mol", regex=True), v / 88.4,
                   np.where(u.str.contains(r"mg/dl|milligram per deciliter", regex=True), v, np.nan))
    lo, hi = K.LABS["CREATININE"]["range_mgdl"]
    return np.where((out >= lo) & (out <= hi), out, np.nan)


def two_tests_apart(lab, flag_col, days=90):
    """P§4: a lab criterion met on 2 tests >= days apart. lab: person_id, date, flag_col. Returns person_id ->
    date the criterion is first met (the second qualifying test)."""
    x = lab.loc[lab[flag_col]].sort_values(["person_id", "date"])
    if x.empty:
        return pd.Series(dtype="datetime64[ns]", index=pd.Index([], name="person_id"))
    first = x.groupby("person_id")["date"].transform("min")
    hit = x.loc[x["date"] >= first + days * DAY]
    return hit.groupby("person_id")["date"].min()


def repeat_codes_apart(ev, days=90):
    """>= 2 code dates >= days apart: person -> date the rule is first met."""
    if ev.empty:
        return pd.Series(dtype="datetime64[ns]", index=pd.Index([], name="person_id"))
    e = ev[["person_id", "date"]].drop_duplicates()
    first = e.groupby("person_id")["date"].transform("min")
    return e.loc[e["date"] >= first + days * DAY].groupby("person_id")["date"].min()


def row_min(cols):
    """Row-wise earliest date of several aligned Series (NaT ignored)."""
    arr = np.column_stack([pd.to_datetime(c).values.astype("datetime64[ns]") for c in cols])
    arr = np.where(np.isnat(arr), np.datetime64("2262-01-01"), arr)
    out = pd.Series(arr.min(axis=1), index=cols[0].index)
    return out.where(out < pd.Timestamp("2262-01-01"))


def near(a, b, lo_days, hi_days, a_col="date", b_col="date"):
    """For each row of a (person_id, a_col): is there a row of b for the same person with
    a_col + lo_days <= b_col <= a_col + hi_days? Returns boolean Series aligned with a."""
    if a.empty or b.empty:
        return pd.Series(False, index=a.index)
    m = a[["person_id", a_col]].reset_index().rename(columns={"index": "_i", a_col: "_a"}).merge(
        b[["person_id", b_col]].rename(columns={b_col: "_b"}).drop_duplicates(), on="person_id")
    k = (m["_b"] >= m["_a"] + lo_days * DAY) & (m["_b"] <= m["_a"] + hi_days * DAY)
    return pd.Series(a.index.isin(m.loc[k, "_i"].unique()), index=a.index)


# =============================================================================
# 4. UTI EPISODES AND RECURRENCE (P§5)
# =============================================================================
def uti_code_events(dx, male=None):
    """dx: person_id, date, set_key (UTI_LOWER / UTI_UPPER / UTI_UNSPEC / PROSTATITIS_ACUTE / UTI_HEMATURIA_FLAG).
    Returns one row per person-date with the code types present that day. Prostatitis counts only when
    D8_include_prostatitis (sensitivity), as upper-equivalent, men only (male: person_id -> bool)."""
    keys = ["UTI_LOWER", "UTI_UPPER", "UTI_UNSPEC"]
    if K.DECISIONS["D8_include_prostatitis"]:
        keys.append("PROSTATITIS_ACUTE")
    e = dx.loc[dx["set_key"].isin(keys + ["UTI_HEMATURIA_FLAG"])].copy()
    if "PROSTATITIS_ACUTE" in keys and male is not None:
        drop = e["set_key"].eq("PROSTATITIS_ACUTE") & ~e["person_id"].map(male).fillna(False).astype(bool)
        e = e.loc[~drop]
    e["upper"] = e["set_key"].isin(["UTI_UPPER", "PROSTATITIS_ACUTE"])
    e["lower"] = e["set_key"].eq("UTI_LOWER")
    e["unspec"] = e["set_key"].eq("UTI_UNSPEC")
    e["hematuria"] = e["set_key"].eq("UTI_HEMATURIA_FLAG")
    e["is_code"] = e["set_key"].isin(keys)
    g = e.groupby(["person_id", "date"], as_index=False)[["upper", "lower", "unspec", "hematuria", "is_code"]].max()
    return g.loc[g["is_code"]].drop(columns="is_code").reset_index(drop=True)


def build_episodes(codes, abx_starts, cultures, gap_days=None, abx_days=None, culture_days=None,
                   require_confirmation=None):
    """P§5.2-5.3. codes: person_id, date, upper, lower, unspec, hematuria (one row per person-date).
    abx_starts: person_id, date (starts of urinary-active systemic antibiotics). cultures: person_id, date.
    Episode = fixed window [start, start + gap_days - 1]; a code >= gap_days after the current start opens a new
    episode. Qualifying = at least one code date with an antibiotic start within +/-abx_days or a culture within
    +/-culture_days (or every episode when require_confirmation is False).
    Returns person_id, ep_start, ep_type (upper/lower), n_code_days, confirmed_abx, confirmed_culture, qualifying,
    hematuria."""
    gap = gap_days or K.DECISIONS["D7_episode_gap_days"]
    ad = K.DECISIONS["D8_confirm_abx_days"] if abx_days is None else abx_days
    cd = K.DECISIONS["D8_confirm_culture_days"] if culture_days is None else culture_days
    req = K.DECISIONS["D8_require_confirmation"] if require_confirmation is None else require_confirmation
    c = codes.sort_values(["person_id", "date"]).reset_index(drop=True)
    if c.empty:
        return pd.DataFrame(columns=["person_id", "ep_start", "ep_type", "n_code_days", "confirmed_abx",
                                     "confirmed_culture", "qualifying", "hematuria"])
    c["c_abx"] = near(c, abx_starts, -ad, ad)
    c["c_cul"] = near(c, cultures, -cd, cd)
    # fixed-window assignment (sequential within person)
    ep_start = np.empty(len(c), dtype="datetime64[ns]")
    pid = c["person_id"].values; dts = c["date"].values.astype("datetime64[ns]")
    gap_ns = np.timedelta64(gap, "D")
    cur_p, cur_s = None, None
    for i in range(len(c)):
        if pid[i] != cur_p or dts[i] >= cur_s + gap_ns:
            cur_p, cur_s = pid[i], dts[i]
        ep_start[i] = cur_s
    c["ep_start"] = ep_start
    E = c.groupby(["person_id", "ep_start"], as_index=False).agg(
        any_upper=("upper", "max"), n_code_days=("date", "nunique"), confirmed_abx=("c_abx", "max"),
        confirmed_culture=("c_cul", "max"), hematuria=("hematuria", "max"))
    E["ep_type"] = np.where(E["any_upper"], "upper", "lower")      # P§5.5 upper dominates; unspecified -> lower (D8b)
    E["qualifying"] = (E["confirmed_abx"] | E["confirmed_culture"]) if req else True
    return E.drop(columns="any_upper")


def index_candidates(E, clean_days=None):
    """P§2.1 / D5: candidate index episodes = qualifying episodes with no episode in the clean_days before them.
    O2: 'any_episode' -> any coded episode blocks; 'qualifying' -> only qualifying episodes block."""
    cd = clean_days or K.DECISIONS["D5_clean_window_days"]
    E = E.sort_values(["person_id", "ep_start"]).copy()
    block = E if K.IMPLEMENTATION["O2_clean_window_counts"] == "any_episode" else E.loc[E["qualifying"]]
    prev = block[["person_id", "ep_start"]].rename(columns={"ep_start": "prev"})
    m = E.loc[E["qualifying"], ["person_id", "ep_start"]].reset_index().merge(prev, on="person_id", how="left")
    m["blocked"] = (m["prev"] < m["ep_start"]) & (m["prev"] >= m["ep_start"] - cd * DAY)
    bl = m.groupby("index")["blocked"].any()
    cand = E.loc[bl.index[~bl.values]].copy()
    return cand.rename(columns={"ep_start": "E1"}).reset_index(drop=True)


def recurrence_in_window(E, idx, window_days=365, within_days=180, k3=3):
    """P§5.4 in the Q1 window [E1, E1 + window_days). idx: person_id, E1 (one row per candidate).
    Returns idx + n_episodes_window, recurrent, single, intermediate, any_upper_window, exposure group."""
    Q = E.loc[E["qualifying"], ["person_id", "ep_start", "ep_type"]]
    m = idx[["person_id", "E1"]].reset_index().merge(Q, on="person_id")
    m = m.loc[(m["ep_start"] >= m["E1"]) & (m["ep_start"] < m["E1"] + window_days * DAY)].sort_values(["index", "ep_start"])
    m["d_prev"] = m.groupby("index")["ep_start"].diff().dt.days
    incl = K.IMPLEMENTATION["O3_recur_180_inclusive"]
    m["close"] = (m["d_prev"] <= within_days) if incl else (m["d_prev"] < within_days)
    agg = m.groupby("index").agg(n_episodes_window=("ep_start", "size"), any_close=("close", "max"),
                                 any_upper_window=("ep_type", lambda s: (s == "upper").any()))
    out = idx.join(agg)
    out["n_episodes_window"] = out["n_episodes_window"].fillna(0).astype(int)
    out["recurrent"] = (out["n_episodes_window"] >= k3) | ((out["n_episodes_window"] >= 2) & out["any_close"].fillna(False).astype(bool))
    out["single"] = out["n_episodes_window"].eq(1)
    out["intermediate"] = out["n_episodes_window"].eq(2) & ~out["recurrent"]
    out["any_upper_window"] = out["any_upper_window"].fillna(False).astype(bool)
    out["q1_group"] = np.select([out["recurrent"], out["single"], out["intermediate"]],
                                ["recurrent", "single", "intermediate"], "none")
    out["uti_type_window"] = np.where(out["any_upper_window"], "any_upper", "lower_only")
    return out.drop(columns="any_close")


# =============================================================================
# 5. ANTIBIOTIC COURSES, INDICATION, PROPHYLAXIS (P§6)
# =============================================================================
def antibiotic_records(drug, gloss):
    """drug rows + glossary -> one row per record and class (a TMP-SMX record is one class B record; a product in two
    classes gives one row per class, flagged). Only counted systemic antibiotics."""
    g = gloss.loc[gloss["abx_counted"], ["concept_id", "abx_classes", "ingredient", "form", "urinary_active",
                                          "prophylaxis_agent", "multi_class"]]
    d = drug.merge(g, left_on="drug_concept_id", right_on="concept_id", how="inner")
    d["abx_class"] = d["abx_classes"].str.split("+")
    d = d.explode("abx_class")
    ds = pd.to_numeric(d["days_supply"], errors="coerce")
    d["rec_end"] = np.where(ds > 0, d["start_date"] + pd.to_timedelta(ds.fillna(0), unit="D"),
                            d["start_date"] + K.DECISIONS["D10_default_course_days"] * DAY)
    d["rec_end"] = pd.to_datetime(d["rec_end"])
    return d.reset_index(drop=True)


def visit_setting(visits):
    v = visits.copy()
    t = v["visit_type"].astype(str).str.lower()
    v["setting"] = np.select([t.str.contains("inpatient"), t.str.contains("emergency")], ["inpatient", "emergency"], "outpatient")
    return v


def build_courses(rec, visits, uti_codes, merge_days=None, ind_days=None):
    """P§6.2, §6.4, §6.5. rec: antibiotic_records output (person_id, abx_class, start_date, rec_end, visit_occurrence_id).
    Inpatient: all records of a class inside one inpatient stay become one course. Then records of the same class
    starting <= merge_days after the previous record's end are merged. Course setting from the first record's visit.
    Indication: urinary when a UTI code (any P§5.1 code) is within +/-ind_days of the course start.
    Returns person_id, abx_class, c_start, c_end, setting, urinary, n_records, with_days_supply."""
    md = merge_days if merge_days is not None else K.DECISIONS["D10_course_merge_days"]
    idd = ind_days if ind_days is not None else K.DECISIONS["D10_indication_days"]
    if rec.empty:
        return pd.DataFrame(columns=["person_id", "abx_class", "c_start", "c_end", "setting", "urinary", "n_records"])
    V = visit_setting(visits)
    r = rec.merge(V[["visit_occurrence_id", "setting", "visit_start_date", "visit_end_date"]], on="visit_occurrence_id", how="left")
    r["setting"] = r["setting"].fillna("outpatient")
    if not K.DECISIONS["D10_include_inpatient"]:
        r = r.loc[r["setting"].ne("inpatient")]
    # one course per class per inpatient stay
    inp = r["setting"].eq("inpatient")
    stay = r.loc[inp].groupby(["person_id", "abx_class", "visit_occurrence_id"], as_index=False).agg(
        start_date=("start_date", "min"), rec_end=("rec_end", "max"), n=("start_date", "size"))
    stay["setting"] = "inpatient"
    out = r.loc[~inp, ["person_id", "abx_class", "start_date", "rec_end", "setting"]].assign(n=1)
    r2 = pd.concat([out, stay.drop(columns="visit_occurrence_id")], ignore_index=True).sort_values(
        ["person_id", "abx_class", "start_date"]).reset_index(drop=True)
    # merge same-class records whose start is <= md days after the running course end
    cid = np.zeros(len(r2), dtype=int)
    p, k, s, e = r2["person_id"].values, r2["abx_class"].values, r2["start_date"].values, r2["rec_end"].values
    cur, run_end, n = None, None, -1
    for i in range(len(r2)):
        key = (p[i], k[i])
        if key != cur or s[i] > run_end + np.timedelta64(md, "D"):
            n += 1; cur = key; run_end = e[i]
        else:
            run_end = max(run_end, e[i])
        cid[i] = n
    r2["cid"] = cid
    Cc = r2.groupby("cid").agg(person_id=("person_id", "first"), abx_class=("abx_class", "first"),
                               c_start=("start_date", "min"), c_end=("rec_end", "max"), setting=("setting", "first"),
                               n_records=("n", "sum")).reset_index(drop=True)
    Cc["urinary"] = near(Cc.rename(columns={"c_start": "date"}), uti_codes, -idd, idd).values
    return Cc


def prophylaxis_periods(rec, uti_codes, longterm_codes, min_days=90, gap_days=14, starts=3, starts_window=120):
    """P§6.5. Prophylaxis per person and agent: (a) >= min_days continuous coverage by one urinary agent allowing
    gaps <= gap_days; (b) >= starts starts of the same agent within starts_window days with no UTI code within +/-3 d
    of any of them; (c) a long-term antibiotic code (Z79.2 / V58.62) with a prescription of the agent within 90 d.
    Returns person_id, agent, p_start, p_end, rule."""
    r = rec.loc[rec["prophylaxis_agent"].notna(), ["person_id", "prophylaxis_agent", "start_date", "rec_end"]].rename(
        columns={"prophylaxis_agent": "agent"}).drop_duplicates().sort_values(["person_id", "agent", "start_date"])
    out = []
    # (a) coverage runs
    for (pid, ag), g in r.groupby(["person_id", "agent"]):
        s0, e0 = None, None
        for s, e in zip(g["start_date"], g["rec_end"]):
            if s0 is None or s > e0 + gap_days * DAY:
                if s0 is not None and (e0 - s0).days >= min_days:
                    out.append((pid, ag, s0, e0, "coverage"))
                s0, e0 = s, e
            else:
                e0 = max(e0, e)
        if s0 is not None and (e0 - s0).days >= min_days:
            out.append((pid, ag, s0, e0, "coverage"))
    # (b) repeated starts without UTI codes
    rs = r.copy()
    rs["no_code"] = ~near(rs.rename(columns={"start_date": "date"}), uti_codes, -3, 3).values
    rs = rs.loc[rs["no_code"]]
    for (pid, ag), g in rs.groupby(["person_id", "agent"]):
        d = g["start_date"].sort_values().tolist()
        for i in range(len(d) - starts + 1):
            if (d[i + starts - 1] - d[i]).days <= starts_window:
                out.append((pid, ag, d[i], d[i + starts - 1] + K.DECISIONS["D10_default_course_days"] * DAY, "repeat_starts"))
                break
    # (c) long-term antibiotic code with a matching prescription
    if longterm_codes is not None and len(longterm_codes):
        m = longterm_codes[["person_id", "date"]].merge(r, on="person_id")
        m = m.loc[(m["start_date"] - m["date"]).abs() <= 90 * DAY]
        for _, x in m.groupby(["person_id", "agent"]).first().reset_index().iterrows():
            out.append((x["person_id"], x["agent"], min(x["date"], x["start_date"]), max(x["date"], x["rec_end"]), "longterm_code"))
    return pd.DataFrame(out, columns=["person_id", "agent", "p_start", "p_end", "rule"])


def apply_prophylaxis_blocks(courses, proph, agent_class=None):
    """P§6.5 course counting: inside a prophylaxis period the agent's courses are replaced by one course per 90-day
    block. agent_class maps agent -> antibiotic class letter (methenamine has none and adds no course)."""
    agent_class = agent_class or {"nitrofurantoin": "A", "trimethoprim": "B", "cephalexin": "D2", "methenamine": None}
    if proph is None or proph.empty:
        return courses.assign(prophylaxis_block=False)
    keep = courses.copy(); keep["prophylaxis_block"] = False
    blocks = []
    for _, p in proph.iterrows():
        cls = agent_class.get(p["agent"])
        if cls is None:
            continue
        inside = keep["person_id"].eq(p["person_id"]) & keep["abx_class"].eq(cls) & \
            keep["c_start"].between(p["p_start"], p["p_end"])
        keep = keep.loc[~inside]
        nb = max(1, int(np.ceil(((p["p_end"] - p["p_start"]).days + 1) / 90)))
        for b in range(nb):
            s = p["p_start"] + b * 90 * DAY
            blocks.append({"person_id": p["person_id"], "abx_class": cls, "c_start": s, "c_end": min(s + 89 * DAY, p["p_end"]),
                           "setting": "outpatient", "urinary": True, "n_records": np.nan, "prophylaxis_block": True})
    out = pd.concat([keep, pd.DataFrame(blocks)], ignore_index=True) if blocks else keep
    return out.drop_duplicates(["person_id", "abx_class", "c_start"]).sort_values(["person_id", "c_start"]).reset_index(drop=True)


# =============================================================================
# 6. KIDNEY CANCER OUTCOME (P§8)
# =============================================================================
def rcc_outcome(dx, px, rcc_drug_starts):
    """P§8.2-8.4. dx: person_id, date, set_key (RCC, RENAL_PELVIS_CA, RENAL_MASS, KIDNEY_DONOR, KIDNEY_TRAUMA,
    STONE_OR_CHRONIC_INF, KIDNEY_TRANSPLANT_STATUS). px: person_id, date, set_key (P§8.3 sets). rcc_drug_starts: person_id, date.
    Returns one row per person with any RCC code: level (confirmed / possible), ambiguous flag, dx_date,
    detection_date, confirm_route."""
    rcc = dx.loc[dx["set_key"].eq("RCC"), ["person_id", "date"]].drop_duplicates()
    if rcc.empty:
        return pd.DataFrame(columns=["person_id", "level", "ambiguous", "dx_date", "detection_date", "confirm_route"])
    # (a) a second RCC code 30-365 days after an earlier RCC code
    two = near(rcc, rcc, 30, 365)
    # (b) a confirming procedure within -90..+180 days of a code, not invalidated
    cpx = px.loc[px["set_key"].isin(K.CONFIRMING_PX), ["person_id", "date"]].drop_duplicates()
    bad_dx = dx.loc[dx["set_key"].isin(["KIDNEY_DONOR", "KIDNEY_TRAUMA", "KIDNEY_TRANSPLANT_STATUS"]), ["person_id", "date"]]
    bad_px = px.loc[px["set_key"].isin(["DONOR_NEPHRECTOMY", "TRANSPLANT_PX"]), ["person_id", "date"]]
    stone = dx.loc[dx["set_key"].eq("STONE_OR_CHRONIC_INF"), ["person_id", "date"]]
    if len(cpx):
        inval = near(cpx, pd.concat([bad_dx, bad_px]), -90, 90)
        inval |= near(cpx, stone, -90, 90) & ~near(cpx, rcc, -90, 90)
        cpx = cpx.loc[~inval]
    proc = near(rcc, cpx, -90, 180)
    # (c) a kidney-cancer drug within 365 days after a code
    drug = near(rcc, rcc_drug_starts, 0, 365)
    rcc = rcc.assign(c_two=two.values, c_proc=proc.values, c_drug=drug.values)
    P = rcc.groupby("person_id").agg(first_code=("date", "min"), n_code_days=("date", "nunique"),
                                     c_two=("c_two", "max"), c_proc=("c_proc", "max"), c_drug=("c_drug", "max"))
    # ambiguous: RCC and renal-pelvis code within 90 days, without a confirming parenchymal procedure
    pel = dx.loc[dx["set_key"].eq("RENAL_PELVIS_CA"), ["person_id", "date"]]
    amb = rcc.loc[near(rcc, pel, -90, 90).values].groupby("person_id").size().gt(0)
    P["ambiguous"] = P.index.map(amb).fillna(False).astype(bool) & ~P["c_proc"]
    P["confirmed"] = (P["c_two"] | P["c_proc"] | P["c_drug"]) & ~P["ambiguous"]
    P["level"] = np.where(P["confirmed"], "confirmed", "possible")
    P["confirm_route"] = np.select([P["c_proc"], P["c_two"], P["c_drug"]], ["procedure", "two_codes", "drug"], "none")
    # dx date = earliest of first RCC code, confirming procedure, biopsy followed by an RCC code within 90 days
    first_px = cpx.groupby("person_id")["date"].min() if len(cpx) else pd.Series(dtype="datetime64[ns]")
    bx = px.loc[px["set_key"].eq("RENAL_BIOPSY"), ["person_id", "date"]].drop_duplicates()
    bx = bx.loc[near(bx, rcc, 0, 90).values].groupby("person_id")["date"].min() if len(bx) else pd.Series(dtype="datetime64[ns]")
    idx_s = P.index.to_series()
    P["dx_date"] = row_min([P["first_code"], first_px.reindex(P.index), bx.reindex(P.index)])
    # detection date: first renal-mass / abnormal-imaging code in the 180 days before dx_date
    rm = dx.loc[dx["set_key"].eq("RENAL_MASS"), ["person_id", "date"]].merge(P["dx_date"].rename("dxd").reset_index(), on="person_id")
    rm = rm.loc[(rm["date"] < rm["dxd"]) & (rm["date"] >= rm["dxd"] - 180 * DAY)]
    P["detection_date"] = rm.groupby("person_id")["date"].min().reindex(P.index)
    return P.reset_index()[["person_id", "level", "ambiguous", "dx_date", "detection_date", "confirm_route", "n_code_days"]]


# =============================================================================
# 7. EHR HISTORY RULE (P§3.2)
# =============================================================================
def history_ok(idx, visit_days, span, t0_col="T0"):
    """idx: person_id, T0 (one row per candidate). visit_days: person_id, date (EHR visit days). span: person_id,
    first_ehr. True when first EHR activity >= 1,095 d before T0 and (D5 rule) >= 1 EHR visit in each of
    [T0-1095, T0-730), [T0-730, T0-365), [T0-365, T0)."""
    f = idx[["person_id", t0_col]].reset_index().merge(span[["person_id", "first_ehr"]], on="person_id", how="left")
    ok = f["first_ehr"] <= f[t0_col] - K.HISTORY_DAYS * DAY
    if K.DECISIONS["D5_history_rule"] == "visit_each_of_3_years":
        m = f[["index", "person_id", t0_col]].merge(visit_days[["person_id", "date"]], on="person_id")
        dd = (m[t0_col] - m["date"]).dt.days
        for lo, hi in ((730, 1095), (365, 730), (0, 365)):     # days before T0: (lo, hi]
            have = m.loc[(dd > lo) & (dd <= hi), "index"].unique()
            ok &= f["index"].isin(have)
    return pd.Series(ok.values, index=idx.index)


# P§3.3 person-level data availability: months with >= 1 EHR condition record / >= 1 EHR drug record (any concept).
ACTIVITY_MONTH_SQL = "WITH " + CAND_CTE + """
SELECT person_id, mon, MAX(has_cond) AS has_cond, MAX(has_drug) AS has_drug FROM (
  SELECT t.person_id, DATE_TRUNC(t.condition_start_date, MONTH) AS mon, 1 AS has_cond, 0 AS has_drug
  FROM `{CDR}.condition_occurrence` t JOIN cand USING (person_id)
  JOIN `{CDR}.condition_occurrence_ext` e ON e.condition_occurrence_id = t.condition_occurrence_id
  WHERE REGEXP_CONTAINS(COALESCE(e.src_id, ''), @ehr_rx)
  UNION ALL
  SELECT t.person_id, DATE_TRUNC(t.drug_exposure_start_date, MONTH), 0, 1
  FROM `{CDR}.drug_exposure` t JOIN cand USING (person_id)
  JOIN `{CDR}.drug_exposure_ext` e ON e.drug_exposure_id = t.drug_exposure_id
  WHERE REGEXP_CONTAINS(COALESCE(e.src_id, ''), @ehr_rx)
) GROUP BY person_id, mon"""


def dmap(keys, s):
    """keys.map(s) for a date Series s (safe when s is empty; pandas cannot cast an empty date map)."""
    s = s[~s.index.duplicated()]
    return pd.Series(pd.to_datetime(s.reindex(keys.values).values), index=keys.index)
