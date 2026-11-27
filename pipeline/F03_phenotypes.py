##### BEGIN F03 phenotypes #####
# UTI episodes (P§5), antibiotic courses and prophylaxis (P§6), RCC levels (P§8), eGFR, ESRD and CKD dates (P§4).
# Descriptive tables here are whole-candidate counts, never split by outcome (blind rule, STUDY-LIFECYCLE stage 7).
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import rutirc as R

C.begin("F03", "phenotypes")
DX, PX, DRUG = C.load_state("dx"), C.load_state("px"), C.load_state("drug")
VIS, LAB, PER = C.load_state("visits"), C.load_state("labs"), C.load_state("persons")
G = C.load_state("drug_glossary", internal=False)
MALE = PER.set_index("person_id")["gender"].astype(str).str.lower().eq("male")


def era(d):
    return np.where(pd.to_datetime(d) < pd.Timestamp(K.ICD10_START), "ICD-9", "ICD-10")


def p_episodes():
    """P§5.2-5.3: 31-day fixed-window episodes; qualifying = code + antibiotic (+/-3 d) or + urine culture (+/-7 d)."""
    codes = R.uti_code_events(DX, MALE)
    rec = R.antibiotic_records(DRUG, G)
    abx = rec.loc[rec["urinary_active"], ["person_id", "start_date"]].rename(columns={"start_date": "date"}).drop_duplicates()
    cul = pd.concat([PX.loc[PX["set_key"].eq("URINE_CULTURE"), ["person_id", "date"]],
                     LAB.loc[LAB["set_key"].eq("URINE_CULTURE_LAB"), ["person_id", "date"]]]).drop_duplicates()
    E = R.build_episodes(codes, abx, cul)
    C.save_state(codes, "uti_codes"); C.save_state(rec, "abx_records"); C.save_state(E, "episodes")
    t = E.assign(era=era(E["ep_start"])).groupby(["era", "ep_type"]).agg(
        n_episodes=("ep_start", "size"), n_qualifying=("qualifying", "sum"), n_by_abx=("confirmed_abx", "sum"),
        n_by_culture=("confirmed_culture", "sum")).reset_index()
    C.show("1", "UTI episodes by coding era and type, and how they were confirmed (P§13 F4, F6)", t,
           counts=["n_episodes", "n_qualifying", "n_by_abx", "n_by_culture"])
    pp = E.loc[E["qualifying"]].groupby("person_id").size()
    C.show("2", "Persons by number of qualifying episodes ever", pp.clip(upper=6).value_counts().sort_index().rename_axis(
        "episodes_ever_capped6").rename("n_persons").reset_index(), counts=["n_persons"], group_cols=[])
C.part("1 UTI episodes", p_episodes)


def p_courses():
    """P§6.2-6.5: courses (7-day merge, one per inpatient stay), indication (+/-3 d), prophylaxis blocks."""
    rec, codes = C.load_state("abx_records"), C.load_state("uti_codes")
    cs = R.build_courses(rec, VIS, codes)
    meth = DRUG.merge(G.loc[G["prophylaxis_only"], ["concept_id"]], left_on="drug_concept_id", right_on="concept_id")
    meth = meth.assign(prophylaxis_agent="methenamine", rec_end=meth["start_date"] + K.DECISIONS["D10_default_course_days"] * R.DAY)
    proph_rec = pd.concat([rec[["person_id", "prophylaxis_agent", "start_date", "rec_end"]],
                           meth[["person_id", "prophylaxis_agent", "start_date", "rec_end"]]], ignore_index=True)
    lt = DX.loc[DX["set_key"].eq("LONGTERM_ABX"), ["person_id", "date"]]
    proph = R.prophylaxis_periods(proph_rec, codes, lt)
    cs2 = R.apply_prophylaxis_blocks(cs, proph)
    C.save_state(cs2, "courses"); C.save_state(proph, "prophylaxis")
    t = cs2.assign(era=era(cs2["c_start"])).groupby(["abx_class", "setting"]).agg(
        n_courses=("c_start", "size"), n_urinary=("urinary", "sum"), n_persons=("person_id", "nunique")).reset_index()
    C.show("3", "Antibiotic courses by class and setting; urinary indication (P§13 F5)", t,
           counts=["n_courses", "n_urinary", "n_persons"])
    C.show("4", "Prophylaxis periods by agent and rule", proph.groupby(["agent", "rule"]).agg(
        n_persons=("person_id", "nunique")).reset_index(), counts=["n_persons"])
    # negative-control exposures (P§10): same course rules, one 'class' per ingredient
    nc = DRUG.merge(G.loc[G["neg_control"], ["concept_id", "ingredient"]], left_on="drug_concept_id", right_on="concept_id")
    nc = nc.assign(abx_class=nc["ingredient"], rec_end=np.where(pd.to_numeric(nc["days_supply"], errors="coerce") > 0,
                   nc["start_date"] + pd.to_timedelta(pd.to_numeric(nc["days_supply"], errors="coerce").fillna(0), unit="D"),
                   nc["start_date"] + K.DECISIONS["D10_default_course_days"] * R.DAY))
    nc["rec_end"] = pd.to_datetime(nc["rec_end"])
    C.save_state(R.build_courses(nc, VIS, codes), "nc_courses")
C.part("2 antibiotic courses", p_courses)


def p_renal():
    """eGFR (CKD-EPI 2021), ESRD date (codes, procedures, eGFR < 15 twice >= 90 d apart), lab CKD date (eGFR < 60)."""
    cr = LAB.loc[LAB["set_key"].eq("CREATININE")].copy()
    cr["scr"] = R.creatinine_mgdl(cr["value_as_number"], cr["unit"])
    bd = R.birth_dates(PER)
    cr["age"] = (cr["date"] - R.dmap(cr["person_id"], bd)).dt.days / 365.25
    cr["female"] = ~cr["person_id"].map(MALE).fillna(False).astype(bool)
    cr = cr.loc[cr["scr"].notna() & cr["age"].ge(18)]
    cr["egfr"] = R.egfr_ckdepi2021(cr["scr"], cr["age"], cr["female"])
    cr["lt15"], cr["lt60"] = cr["egfr"] < 15, cr["egfr"] < 60
    C.save_state(cr[["person_id", "date", "egfr"]], "egfr")
    esrd = pd.concat([DX.loc[DX["set_key"].eq("ESRD_DX"), ["person_id", "date"]],
                      PX.loc[PX["set_key"].isin(["DIALYSIS_PX", "TRANSPLANT_PX"]), ["person_id", "date"]],
                      R.two_tests_apart(cr, "lt15").rename("date").reset_index()])
    C.save_state(esrd.groupby("person_id")["date"].min().rename("esrd_date").reset_index(), "esrd")
    C.save_state(R.two_tests_apart(cr, "lt60").rename("ckd_lab_date").reset_index(), "ckd_lab")
    t = pd.DataFrame([{"n_creatinine_tests": len(cr), "n_persons": cr["person_id"].nunique(),
                       "share_unit_unusable": float(np.mean(np.isnan(R.creatinine_mgdl(
                           LAB.loc[LAB["set_key"].eq("CREATININE"), "value_as_number"], LAB.loc[LAB["set_key"].eq("CREATININE"), "unit"]))))}])
    C.show("5", "Creatinine tests usable for eGFR", t, counts=["n_creatinine_tests", "n_persons"])
C.part("3 eGFR, ESRD, CKD", p_renal)


def p_outcome():
    """P§8: RCC confirmed / possible / ambiguous, diagnosis and detection dates. Whole-candidate totals only (F1)."""
    rd = DRUG.merge(G.loc[G["rcc_drug"], ["concept_id"]], left_on="drug_concept_id", right_on="concept_id")[["person_id", "start_date"]]
    O = R.rcc_outcome(DX, PX, rd.rename(columns={"start_date": "date"}))
    C.save_state(O, "rcc")
    t = O.groupby(["level", "confirm_route"]).size().rename("n_persons").reset_index()
    C.show("6", "Persons with any RCC code, by level and confirmation route (P§13 F1; candidates, all dates)", t,
           counts=["n_persons"], group_cols=[])
    t2 = O.assign(era=era(O["dx_date"])).groupby(["era", "level"]).size().rename("n_persons").reset_index()
    C.show("7", "RCC persons by coding era of diagnosis date (P§13 F6)", t2, counts=["n_persons"], group_cols=[])
C.part("4 kidney cancer outcome", p_outcome)
C.end()
##### END F03 #####
