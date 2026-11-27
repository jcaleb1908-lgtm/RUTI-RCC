##### BEGIN F02 extraction #####
# Candidate persons = anyone with a P§5.1 UTI code or an RCC code. Pulls only the concept sets resolved in F00.
# Everything saved here is participant-level: _INTERNAL, stays in the run folder, never copied out.
import pandas as pd
import aou_core as C
import _config as K
import rutirc as R

C.begin("F02", "extraction")
S = C.load_state("code_sets", internal=False)
G = C.load_state("drug_glossary", internal=False)
ATC = C.load_state("atc_drugs", internal=False)
EHR_RX = K.IMPLEMENTATION["O5_ehr_src_regex"]
CAND_KEYS = ["UTI_LOWER", "UTI_UPPER", "UTI_UNSPEC", "PROSTATITIS_ACUTE", "RCC"]
cand_ids = S.loc[S["set_key"].isin(CAND_KEYS) & S["kind"].eq("DX"), "concept_id"].astype(int).unique().tolist()
CAND = C.A("cand_ids", "INT64", cand_ids)
MG = K.RUN["max_gib"]


def p_dx():
    """All diagnosis code sets (P§4, 5, 8, 9, 10) from condition, observation, procedure and measurement tables."""
    dx_sets = S.loc[S["kind"].eq("DX"), ["concept_id", "set_key", "era"]]
    d = C.q(R.DX_EVENTS_SQL, [CAND, C.A("ids", "INT64", dx_sets["concept_id"].astype(int).unique().tolist())], "dx events", max_gib=MG)
    d["date"] = C.to_date(d["event_date"])
    d = d.merge(dx_sets, left_on="source_concept_id", right_on="concept_id").drop(columns=["concept_id", "event_date"])
    C.save_state(d, "dx")
    C.show("1", "Diagnosis records by table of origin (all sets)", d.groupby("tbl").size().rename("n_records").reset_index())
C.part("1 diagnoses", p_dx)


def p_px():
    """Procedure code sets (P§4, 5.3, 8.3, 9.4): source and standard procedure concept, plus measurement/observation."""
    px_sets = S.loc[S["kind"].eq("PX"), ["concept_id", "set_key"]]
    d = C.q(R.PX_EVENTS_SQL, [CAND, C.A("ids", "INT64", px_sets["concept_id"].astype(int).unique().tolist())], "px events", max_gib=MG)
    d["date"] = C.to_date(d["event_date"])
    d = d.merge(px_sets, on="concept_id").drop(columns="event_date").drop_duplicates(["person_id", "date", "set_key"])
    C.save_state(d, "px")
C.part("2 procedures", p_px)


def p_drug():
    """Study antibiotics, negative-control drugs, kidney-cancer drugs, methenamine and ATC covariate groups."""
    want = G.loc[G["abx_counted"] | G["neg_control"] | G["rcc_drug"] | G["prophylaxis_only"], "concept_id"]
    ids = sorted(set(want.astype(int)) | set(ATC["concept_id"].astype(int)))
    d = C.q(R.DRUG_SQL, [CAND, C.A("ids", "INT64", ids)], "drug exposures", max_gib=MG)
    d["start_date"] = C.to_date(d["start_date"]); d["end_date"] = C.to_date(d["end_date"])
    C.save_state(d, "drug")
C.part("3 drugs", p_drug)


def p_visits():
    """EHR visits (P§3.2) with type and src_id; per-person EHR span from all EHR-sourced tables."""
    v = C.q(R.VISITS_SQL, [CAND], "visits", max_gib=MG)
    for c in ("visit_start_date", "visit_end_date"):
        v[c] = C.to_date(v[c])
    v["is_ehr"] = v["src_id"].fillna("").str.contains(EHR_RX, regex=True)
    C.save_state(v, "visits")
    sp = C.q(R.EHR_SPAN_SQL, [CAND, C.S("ehr_rx", "STRING", EHR_RX)], "EHR span", max_gib=MG)
    sp["first_ehr"] = C.to_date(sp["first_ehr"]); sp["last_ehr"] = C.to_date(sp["last_ehr"])
    C.save_state(sp, "ehr_span")
    am = C.q(R.ACTIVITY_MONTH_SQL, [CAND, C.S("ehr_rx", "STRING", EHR_RX)], "EHR activity months", max_gib=MG)
    am["mon"] = C.to_date(am["mon"])
    C.save_state(am, "activity_months")
C.part("4 visits and EHR span", p_visits)


def p_labs():
    """Creatinine, HbA1c, BMI, urine culture (LOINC sets from F00)."""
    lab_sets = S.loc[S["kind"].eq("LAB"), ["concept_id", "set_key"]]
    d = C.q(R.LABS_SQL, [CAND, C.A("ids", "INT64", lab_sets["concept_id"].astype(int).tolist())], "labs", max_gib=MG)
    d["date"] = C.to_date(d["measurement_date"])
    d = d.merge(lab_sets, left_on="measurement_concept_id", right_on="concept_id").drop(columns=["concept_id", "measurement_date"])
    C.save_state(d, "labs")
C.part("5 labs", p_labs)


def p_person():
    """Person, death (primary record), enrollment (primary consent), survey items."""
    p = C.q(R.PERSON_SQL, [CAND], "persons")
    C.save_state(p, "persons")
    C.save_state(C.death_dates().rename("death_date").reset_index(), "death")
    C.save_state(C.enrollment_dates().rename("enroll_date").reset_index(), "enroll")
    inv = C.load_state("survey_inventory", internal=False, required=False)
    qcodes = list(K.SURVEY_ITEMS)
    if inv is not None:
        extra = inv.loc[inv["question"].str.contains(K.HISPANIC_SUBGROUP_QUESTION_RX + "|" + K.BIRTHPLACE_QUESTION_RX,
                                                     case=False, regex=True), "question_code"].unique().tolist()
        qcodes = sorted(set(qcodes) | set(extra))
    s = C.q(R.SURVEY_SQL, [CAND, C.A("qcodes", "STRING", qcodes)], "survey answers", max_gib=MG)
    s["date"] = C.to_date(s["observation_date"])
    C.save_state(s, "survey")
    C.show("2", "Candidate persons", pd.DataFrame([{"n": p["person_id"].nunique()}]))
C.part("6 persons, death, enrollment, survey", p_person)
C.end()
##### END F02 #####
