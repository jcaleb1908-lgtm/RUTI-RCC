##### BEGIN F01 data inventory and quality (P§13 F7, §16.2) #####
# Whole-CDR aggregates only. Answers the open questions the design depends on (M6, M7) before any cohort is built.
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import rutirc as R

C.begin("F01", "data inventory and quality")
EHR_RX = K.IMPLEMENTATION["O5_ehr_src_regex"]


def p_dob():
    """16.2 Q1: is the full date of birth released?"""
    d = C.q(R.DOB_SQL, [], "date of birth precision")
    C.show("1", "Date-of-birth precision (all persons)", d, counts=["persons", "n_day_missing", "n_day_15", "n_birth_datetime_missing"])
C.part("1 date of birth", p_dob)


def p_src():
    """16.2 Q2: is a data-partner identifier (src_id) available? Patterns only; sites are never named (AL7)."""
    d = C.q(R.SRC_PATTERN_SQL, [], "src_id patterns")
    d["counted_as_ehr"] = d["src_pattern"].str.contains(EHR_RX, regex=True)
    C.show("2", "src_id patterns on visit records (digits replaced by #)", d, counts=["persons"])
C.part("2 data partner identifier", p_src)


def p_drugs():
    """16.2 Q3: days supply completeness; prescription vs inpatient administration (study antibiotics only)."""
    G = C.load_state("drug_glossary", internal=False)
    ids = G.loc[G["abx_counted"], "concept_id"].astype(int).tolist()
    d = C.q(R.DRUG_QUALITY_SQL, [C.A("ids", "INT64", ids)], "antibiotic record quality", max_gib=K.RUN["max_gib"])
    C.show("3", "Study antibiotic records by drug type: completeness of days supply, end date, route, visit link", d,
           counts=["n_records", "persons"])
C.part("3 drug record quality", p_drugs)


def p_culture():
    """16.2 Q4: urine culture measurements located by name; are results structured (value_as_concept)?"""
    d = C.q(R.LAB_NAME_SEARCH_SQL, [C.S("rx", "STRING", K.LABS["URINE_CULTURE_LAB"]["name_search"])], "urine culture by name",
            max_gib=K.RUN["max_gib"])
    C.show("4", "Measurement concepts named like urine culture", d.sort_values("persons", ascending=False).head(60),
           counts=["persons"])
C.part("4 urine culture", p_culture)


def p_death():
    """16.2 Q5: death record sources."""
    d = C.q(R.DEATH_SOURCE_SQL, [], "death sources")
    C.show("5", "aou_death rows by source pattern and primary flag", d, counts=["persons"])
C.part("5 death source", p_death)


def p_survey():
    """16.2 Q7: Basics / Lifestyle items and answers by name (fills CROSSWALKS and checks SURVEY_ITEMS codes)."""
    d = C.q(R.SURVEY_INVENTORY_SQL, [C.S("rx", "STRING", K.SURVEY_SEARCH_RX)], "survey inventory", max_gib=K.RUN["max_gib"])
    C.save_state(d.drop(columns="persons"), "survey_inventory", internal=False)
    C.show("6", "Survey questions and answers matching race, ethnicity, birthplace, sex, education, income, insurance, smoking, ZIP",
           d.sort_values(["question_code", "answer_code"]), counts=["persons"], group_cols=["question_code"], max_rows=600)
    missing = [c for c in K.SURVEY_ITEMS if c not in set(d["question_code"])]
    C.show("7", "SURVEY_ITEMS codes not found in this CDR (fix _config.py before F02)", pd.DataFrame({"question_code": missing}),
           counts=[])
    for var, cw in K.CROSSWALKS.items():
        q = [k for k, v in K.SURVEY_ITEMS.items() if v == var]
        ans = set(d.loc[d["question_code"].isin(q), "answer_code"].dropna())
        C.show(f"8_{var}", f"{var}: answers in the CDR missing from CROSSWALKS (would be 'Not mapped')",
               pd.DataFrame({"answer_code": sorted(ans - set(cw))}), counts=[])
    sub = d.loc[d["question"].str.contains(K.HISPANIC_SUBGROUP_QUESTION_RX, case=False, regex=True)
                | d["answer"].fillna("").str.contains(r"puerto", case=False, regex=True)]
    C.show("9", "Hispanic subgroup / Puerto Rico answers (P§11.1; is the Puerto Rican subgroup released?)", sub,
           counts=["persons"], group_cols=["question_code"])
C.part("6 survey inventory", p_survey)


def p_zip():
    """P§11.1, 16.2 Q6: participants with a Puerto Rico ZIP3."""
    code = next((k for k, v in K.SURVEY_ITEMS.items() if v == "zip3"), None)
    d = C.q(R.ZIP3_SQL, [C.S("zip_code", "STRING", code), C.A("zips", "STRING", K.PR_ZIP3)], "PR zip3")
    C.show("10", "Participants by Puerto Rico ZIP3 prefix", d, counts=["persons"], group_cols=[])
C.part("7 Puerto Rico ZIP3", p_zip)


def p_site_year():
    """P§3.3: drug records per 1,000 visits by site and year. The per-site table stays INTERNAL; only its
    distribution is shown."""
    d = C.q(R.SITE_YEAR_SQL, [C.S("ehr_rx", "STRING", EHR_RX)], "site-year drug completeness", max_gib=K.RUN["max_gib"])
    d["drugs_per_1000_visits"] = 1000 * d["n_drugs"] / d["n_visits"].clip(lower=1)
    C.save_state(d, "site_year")
    p10 = float(np.percentile(d["drugs_per_1000_visits"], 10)) if len(d) else np.nan
    t = pd.DataFrame([{"site_years": len(d), "sites": d["src_id"].nunique(),
                       "ratio_p10_p25_p50_p75": C.pctl(d["drugs_per_1000_visits"], (10, 25, 50, 75)),
                       "site_years_below_p10": int((d["drugs_per_1000_visits"] < p10).sum()),
                       "site_years_no_drugs": int(d["n_drugs"].eq(0).sum())}])
    C.show("11", "Site-year drug completeness (distribution only; sites not shown)", t, counts=[])
C.part("8 site-year completeness", p_site_year)


def p_obs_period():
    """16.2 Q8: observation_period exists (not used for the history rule until its construction is known)."""
    C.show("12", "observation_period size", C.q(R.OBS_PERIOD_SQL, [], "observation period"), counts=["persons"])
C.part("9 observation period", p_obs_period)
C.end()
##### END F01 #####
