"""_config.py - RUTI-RCC single source of truth (ANALYSIS-STANDARD.md section 3).

Study: recurrent UTI, its antibiotics, and renal cell carcinoma (RCC) in All of Us (Controlled Tier, CDR v9).
Protocol: docs/PROTOCOL-ruti-rcc-v0.1.md (Pre-cohort protocol v0.1, 19 Nov 2026, draft for PI review).

Rules for this file
  * Every definition comes from a protocol section, cited in the comment next to it (P§x.y).
  * Every open PI decision (D1-D21) is a named switch with the protocol's PROPOSED value as default.
  * No All of Us concept IDs are typed here (protocol "How to read"; MISTAKES M2). Code sets are written as
    source codes (ICD-10-CM, ICD-9-CM, CPT4, ICD-10-PCS, ICD-9 procedures, LOINC, ATC) or RxNorm ingredient
    names. F00 resolves them inside the Workbench and prints every concept it reached, by name.
  * Code syntax: "N30.0*" = the code and every code under it; "N10" = that exact code.
  * Bump PIPELINE_VERSION on any definition change (printed in every table footnote).
"""

PROJECT = "RUTI-RCC"
PIPELINE_VERSION = "RUTI-RCC-2026.11.28-v0.2"      # protocol v0.1 decisions, PI-approved 2026-11-28
PROTOCOL_VERSION = "Pre-cohort protocol v0.1 (2026-11-19, revised 2026-11-28); D1-D21 PI-approved 2026-11-28"
CDR_EXPECTED = "C2025Q4R6"                         # P§14: v9 Controlled Tier; data through 2025-01-01
CDR_CUTOFF = "2025-01-01"                          # M17: every date checked against it

# ---------------------------------------------------------------------------------------------------------
# PI DECISIONS (P§16). The protocol's "Proposed" values, approved by the PI on 2026-11-28. Change only by AMENDMENT.
# ---------------------------------------------------------------------------------------------------------
DECISIONS = {
    "D2_q1_time_zero": "landmark_365",          # alt: "sequential_trials", "time_varying"
    "D3_q2_reference": "A_nitrofuran_fosfomycin",
    "D4_t0_before_enrollment": True,            # primary allows T0 before enrollment; delayed entry = sensitivity
    "D5_history_rule": "visit_each_of_3_years", # alt: "record_span"
    "D5_t0_window": ("2008-01-01", "2022-12-31"),
    "D5_clean_window_days": 365,                # no UTI episode in the 365 days before E1
    "D6_hereditary_before_t0_only": True,
    "D7_episode_gap_days": 31,                  # sensitivity 14
    "D8_confirm_abx_days": 3,                   # code + urinary-active antibiotic within +/-3 d
    "D8_confirm_culture_days": 7,               # or code + urine culture within +/-7 d
    "D8_require_confirmation": True,            # False = code-only sensitivity
    "D8_n39_counts_as_lower": True,
    "D8_include_prostatitis": False,            # True = sensitivity: acute prostatitis as upper-equivalent in men
    "D9_primary_exposure": "recurrent_any_type",
    "D10_course_merge_days": 7,
    "D10_default_course_days": 7,               # when days supply is missing
    "D10_indication_days": 3,                   # sensitivity 7
    "D10_baseline_abx_years": 3,
    "D10_include_inpatient": True,
    "D11_lag_years": 2,
    "D11_latency_years": 2,
    "D12_outcome_level": "confirmed",           # alt: "confirmed_or_possible"
    "D13_both_sexes": True,
    "D14_nco": ["NCO_ADRENAL_BENIGN", "NCO_CATARACT"],
    "D15_subgroup_min_events": 100,
    "D15_subgroup_min_per_cell": 21,
    "D16_missing": "mice",                      # m = 20 (P§12.2 #8)
    "D17_go_events": 300,
    "D17_reduced_events": 100,
    "D20_censor_esrd": True,
}

# Implementation choices the protocol leaves open (README "Open items" O#). Defaults chosen here, flagged for the PI.
IMPLEMENTATION = {
    "O1_index_choice": "first_eligible",        # E1 = earliest clean qualifying episode whose T0 passes eligibility;
                                                # alt "first_ever" = earliest clean qualifying episode only
    "O2_clean_window_counts": "any_episode",    # a coded (even unconfirmed) episode in the prior 365 d blocks E1
    "O3_recur_180_inclusive": True,             # ">=2 episodes within 180 days" read as start difference <= 180 d
    "O4_q2_multi_class_scope": "A-G",           # episodes starting >=2 of classes A-G within +/-3 d are excluded
    "O5_ehr_src_regex": r"^EHR site",           # src_id values counted as EHR (checked in F01)
    "O6_dob_day_if_missing": 15,                # P§3.1
}

AGE_MIN = 40                                    # P§3.1, completed years at T0
HISTORY_DAYS = 1095                             # P§3.2
WASHOUT_Q2_DAYS = 90                            # P§6.6
FOLLOWUP_CUTOFF = CDR_CUTOFF
ICD10_START = "2015-10-01"                      # P§3.4 coding era
ICD10_ERA_ONLY_T0 = "2018-10-01"                # P§12.2 #15

# ---------------------------------------------------------------------------------------------------------
# DIAGNOSIS CODE SETS (ICD-10-CM / ICD-9-CM source codes), P§4, §5, §8, §9, Appendix A.1
# role: E = exposure/episode, O = outcome, X = exclusion, C = covariate, N = negative control, D = detection
# ---------------------------------------------------------------------------------------------------------
DX = {
    # P§5.1 UTI episode codes
    "UTI_LOWER": dict(role="E", icd10=["N30.00", "N30.01", "N30.80", "N30.81", "N30.90", "N30.91", "O23.1*", "O86.22"],
                      icd9=["595.0", "595.89", "595.9"], rx=r"cystitis|bladder|urinary"),
    "UTI_UPPER": dict(role="E", icd10=["N10", "N12", "N13.6", "N15.1", "O23.0*", "O86.21"],
                      icd9=["590.10", "590.11", "590.2", "590.80", "590.9"],
                      rx=r"pyelonephritis|nephritis|pyonephrosis|abscess|kidney|infection"),
    "UTI_UNSPEC": dict(role="E", icd10=["N39.0", "O23.4*", "O86.20", "O86.29"], icd9=["599.0"],
                       rx=r"urinary tract|urinary system|genitourinary|infection"),
    "PROSTATITIS_ACUTE": dict(role="E", icd10=["N41.0"], icd9=["601.0"], rx=r"prostatitis"),
    "UTI_HEMATURIA_FLAG": dict(role="C", icd10=["N30.01", "N30.81", "N30.91"], icd9=[], rx=r"hematuria"),
    "HX_UTI": dict(role="C", icd10=["Z87.440"], icd9=["V13.02"], rx=r"urinary"),
    # P§8.1 outcome
    "RCC": dict(role="O", icd10=["C64*"], icd9=["189.0"], rx=r"kidney|renal"),
    "RENAL_PELVIS_CA": dict(role="O", icd10=["C65*"], icd9=["189.1"], rx=r"pelvis"),
    "URETER_CA": dict(role="O", icd10=["C66*"], icd9=["189.2"], rx=r"ureter"),
    "BLADDER_CA": dict(role="O", icd10=["C67*"], icd9=["188*"], rx=r"bladder|urachus|ureteric orifice"),
    "RENAL_MASS": dict(role="D", icd10=["D30.0*", "D41.0*", "R93.42*"], icd9=[], rx=r"kidney|renal"),
    # P§4 exclusions
    "HX_KIDNEY_CA": dict(role="X", icd10=["Z85.52*"], icd9=["V10.52", "V10.53"], rx=r"kidney|renal"),
    "OTHER_CANCER": dict(role="X", icd10=[f"C{i:02d}*" for i in range(0, 97)] + ["C4A*", "C7A*", "C7B*", "D09.0", "D09.1*", "Z85*"],
                         icd10_exclude=["C44*", "Z85.828"],
                         icd9=[f"{i}*" for i in range(140, 209)] + ["209.0*", "209.1*", "209.2*", "209.3*", "V10*"],
                         icd9_exclude=["173*", "V10.83"], rx=r"malignan|neoplasm|cancer|carcinoma|lymphoma|leukemia|myeloma|tumor|sarcoma|melanoma|history|"
                            r"hodgkin|mesothelioma|plasmacytoma|myelofibrosis|erythremia|heavy chain|mycosis fungoides|"
                            r"sezary|histiocytosis|letterer|reticuloendotheliosis|macroglobulinemia|mastocytosis|"
                            r"immunoproliferative|myelodysplastic|hepatoblastoma|cd30-positive|panmyelosis"),
    "ESRD_DX": dict(role="X", icd10=["N18.6", "Z99.2", "Z49*", "Z91.15*", "Z94.0", "T86.1*"],
                    icd9=["585.6", "V45.11", "V45.12", "V56*", "V42.0", "996.81"],
                    rx=r"end stage|dialysis|transplant|dependence on renal|noncompliance"),
    "VHL": dict(role="X", icd10=["Q85.83", "Q85.8"], icd9=["759.6"], rx=r"hippel|phakomatos|hamartos"),
    "TSC": dict(role="X", icd10=["Q85.1"], icd9=["759.5"], rx=r"tuberous"),
    "PKD": dict(role="X", icd10=["Q61.11", "Q61.19", "Q61.2", "Q61.3"], icd9=["753.12", "753.13", "753.14"], rx=r"polycystic"),
    "CATHETER_ANY": dict(role="X", icd10=["Z93.5*"], icd9=["V44.5*"], rx=r"cystostomy|vesicostomy"),          # any one code
    "CATHETER_REPEAT": dict(role="X", icd10=["Z96.0", "Z46.6", "T83.51*"], icd9=["V53.6", "996.64"],
                            rx=r"urogenital|urinary|catheter"),                                  # >=2 codes >=90 d apart
    "NEUROGENIC_BLADDER": dict(role="X", icd10=["N31*"], icd9=["596.54", "344.61"], rx=r"neurogenic|neuromuscular|bladder"),
    "CONGENITAL_URO": dict(role="X", icd10=["Q60*", "Q62*", "Q63*", "Q64*", "Q61.4", "Q61.5", "Q61.8", "Q61.9"],
                           icd9=["753.0", "753.15", "753.16", "753.17", "753.19", "753.2*", "753.3", "753.4", "753.5",
                                 "753.6", "753.7", "753.8", "753.9"], rx=r"kidney|renal|ureter|bladder|urinary|congenital|cystic|anomal|agenesis|hypoplasia|dysplasia|exstrophy|urachus|cloaca|duplication|obstructive"),
    # P§4.1, §9.4 covariates
    "STONES": dict(role="C", icd10=["N20*", "Z87.442"], icd9=["592*", "V13.01"], rx=r"calcul|stone|urolith"),
    "HIV": dict(role="C", icd10=["B20", "Z21"], icd9=["042"], rx=r"immunodeficiency|HIV"),
    "IMMUNO_DX": dict(role="C", icd10=["D84.821", "Z79.52", "Z79.6*"], icd9=["V58.65", "V87.46", "279.3"],
                      rx=r"immun|steroid|immunosuppress|long term \(current\) use"),
    "DIABETES": dict(role="C", icd10=["E08*", "E09*", "E10*", "E11*", "E13*"], icd9=["250*"], rx=r"diabet"),
    "HYPERTENSION": dict(role="C", icd10=["I10", "I11*", "I12*", "I13*", "I15*", "I16*"], icd9=["401*", "402*", "403*", "404*", "405*"],
                         rx=r"hypertens"),
    "CKD": dict(role="C", icd10=["N18.1", "N18.2", "N18.30", "N18.31", "N18.32", "N18.4", "N18.5", "N18.9"],
                icd9=["585.1", "585.2", "585.3", "585.4", "585.5", "585.9"], rx=r"chronic kidney"),
    "SMOKING_DX": dict(role="C", icd10=["F17.2*", "Z72.0", "Z87.891"], icd9=["305.1", "V15.82"], rx=r"nicotine|tobacco|smok"),
    "OBESITY_DX": dict(role="C", icd10=["E66*", "Z68.3*", "Z68.4*"], icd9=["278.00", "278.01", "V85.3*", "V85.4*"],
                       rx=r"obes|body mass"),
    "HEMATURIA": dict(role="C", icd10=["R31*", "N02*"], icd9=["599.70", "599.71", "599.72"], rx=r"hematuria"),
    "VUR": dict(role="C", icd10=["N13.70", "N13.71", "N13.72*", "N13.73*", "Q62.7"], icd9=["593.7*"], rx=r"reflux"),
    "CHRONIC_PYELO": dict(role="C", icd10=["N11*"], icd9=["590.00", "590.01"], rx=r"pyelonephritis|tubulo"),
    "LONGTERM_ABX": dict(role="C", icd10=["Z79.2"], icd9=["V58.62"], rx=r"antibiotic"),
    "PCN_ALLERGY": dict(role="C", icd10=["Z88.0"], icd9=[], rx=r"penicillin"),
    "RESISTANT_ORG": dict(role="C", icd10=["Z16*"], icd9=["V09*"], rx=r"resist"),   # P§9.1 ResistantOrg0 [proposal: codes not in protocol]
    # P§8.3 procedure-confirmation exclusions
    "KIDNEY_DONOR": dict(role="C", icd10=["Z52.4"], icd9=["V59.4"], rx=r"donor"),
    "KIDNEY_TRANSPLANT_STATUS": dict(role="C", icd10=["Z94.0", "T86.1*"], icd9=["V42.0", "996.81"], rx=r"transplant"),
    "KIDNEY_TRAUMA": dict(role="C", icd10=["S37.0*"], icd9=["866*"], rx=r"kidney|renal"),
    "STONE_OR_CHRONIC_INF": dict(role="C", icd10=["N20*", "N11*", "N13.6", "N15.1"], icd9=["592*", "590.0*", "590.2"],
                                 rx=r"calcul|stone|pyelonephritis|pyonephrosis|abscess|tubulo"),
    # P§10, D14 negative-control outcomes (+ alternatives, computed but not primary)
    "NCO_ADRENAL_BENIGN": dict(role="N", icd10=["D35.00", "D35.01", "D35.02"], icd9=["227.0"], rx=r"adrenal"),
    # Only codes whose ICD title is a cataract: 366.11 "Pseudoexfoliation of lens capsule" is not (decision 2026-11-28).
    "NCO_CATARACT": dict(role="N", icd10=["H25*"], icd9=["366.1*"], icd9_exclude=["366.11"], rx=r"cataract|nuclear sclerosis"),
    "NCO_MELANOMA": dict(role="N", icd10=["C43*"], icd9=["172*"], rx=r"melanoma"),
}

# ---------------------------------------------------------------------------------------------------------
# PROCEDURES (CPT4 / ICD-10-PCS / ICD-9 procedure source codes), P§4, §5.3, §8.3, §9.4, Appendix A.2
# ---------------------------------------------------------------------------------------------------------
PX = {
    "NEPHRECTOMY": dict(cpt=["50220", "50225", "50230", "50545", "50546"],
                        pcs=["0TT00ZZ", "0TT04ZZ", "0TT04ZG", "0TT10ZZ", "0TT14ZZ", "0TT14ZG", "0TT20ZZ", "0TT24ZZ", "0TT24ZG"],
                        icd9p=["55.52", "55.54"], rx=r"nephrectomy|kidney|resection"),
    # Nephroureterectomy (kidney + whole ureter; mostly for renal pelvis / ureter cancer). Confirms RCC in the primary
    # analysis; a sensitivity run drops it (decision 2026-11-28). ICD-10-PCS has no single code for it.
    "NEPHROURETERECTOMY": dict(cpt=["50234", "50236", "50548"], pcs=[], icd9p=["55.51"],
                               rx=r"nephrectomy|nephroureterectomy|ureterectomy"),
    "PARTIAL_NEPHRECTOMY": dict(cpt=["50240", "50543"], pcs=["0TB00ZZ", "0TB03ZZ", "0TB04ZZ", "0TB10ZZ", "0TB13ZZ", "0TB14ZZ"],
                                icd9p=["55.4"], rx=r"nephrectomy|kidney|excision"),
    "RENAL_ABLATION": dict(cpt=["50250", "50542", "50592", "50593"], pcs=["0T500ZZ", "0T503ZZ", "0T504ZZ", "0T510ZZ", "0T513ZZ", "0T514ZZ"],
                           icd9p=["55.32", "55.33", "55.34", "55.35", "55.39"], rx=r"ablation|destruction|kidney|renal|cryo"),
    "RENAL_BIOPSY": dict(cpt=["50200", "50205"], pcs=["0TB03ZX", "0TB13ZX"], icd9p=["55.23", "55.24"], rx=r"biopsy|kidney|excision"),
    "DONOR_NEPHRECTOMY": dict(cpt=["50300", "50320", "50547"], pcs=[], icd9p=[], rx=r"donor"),
    "TRANSPLANT_PX": dict(cpt=["50340", "50360", "50365"], pcs=["0TY00Z0", "0TY00Z1", "0TY00Z2", "0TY10Z0", "0TY10Z1", "0TY10Z2"],
                          icd9p=["55.69"], rx=r"transplant|nephrectomy|kidney"),
    "DIALYSIS_PX": dict(cpt=["90935", "90937", "90945", "90947"], pcs=[], icd9p=["39.95", "54.98"], rx=r"dialysis|hemodialysis|hemofiltration"),
    "URINE_CULTURE": dict(cpt=["87086", "87088"], pcs=[], icd9p=[], rx=r"culture|urine"),
    "IMAGING_ABD": dict(cpt=["74150", "74160", "74170", "74176", "74177", "74178", "74181", "74182", "74183",
                             "76700", "76705", "76770", "76775", "72192", "72193", "72194"],
                        pcs=[], icd9p=["87.71", "88.01"], rx=r"abdom|pelvi|retroperiton|tomograph|ultrasound|magnetic|CT|MRI|echography|kidney"),
    "STONE_PX": dict(cpt=["50080", "50081", "50590", "52352", "52353", "52356"], pcs=[], icd9p=[], rx=r"lithotr|nephrostolith|calcul|ureteroscop"),
    "CYSTOSCOPY": dict(cpt=["52000"], pcs=[], icd9p=[], rx=r"cystourethroscopy|cystoscop"),
}
# Kidney-cancer-confirming procedures (P§8.3: biopsy supports the date but never confirms alone)
CONFIRMING_PX = ["NEPHRECTOMY", "NEPHROURETERECTOMY", "PARTIAL_NEPHRECTOMY", "RENAL_ABLATION"]
CONFIRMING_PX_SENS = [k for k in CONFIRMING_PX if k != "NEPHROURETERECTOMY"]   # sensitivity: without nephroureterectomy

# ---------------------------------------------------------------------------------------------------------
# LABS (LOINC codes; name-checked in F00), P§4, §4.1, §9.4
# The LOINC codes below are proposals to confirm by name in F00 (the protocol names the labs, not codes).
# ---------------------------------------------------------------------------------------------------------
LABS = {
    "CREATININE": dict(loinc=["2160-0", "38483-4"], rx=r"creatinine", range_mgdl=(0.1, 20.0)),
    "HBA1C": dict(loinc=["4548-4", "17856-6"], rx=r"a1c|hemoglobin a1", range=(3.0, 20.0)),
    "BMI": dict(loinc=["39156-5"], rx=r"body mass", range=(12.0, 90.0)),
    # P§5.3 "a LOINC urine-culture measurement found by name": from the F01.4 inventory (2026-11-28 run). Reflex-culture
    # panels (58077-9, 57019-2, LP402549-2) and virus / mycobacteria cultures are left out. SNOMED codes = O9 proposal.
    "URINE_CULTURE_LAB": dict(loinc=["630-4", "13315-7", "17970-5", "17971-3", "17972-1", "17973-9", "17974-7",
                                     "44847-2", "44849-8"],
                              snomed=["104193001", "104194007", "117010004"],
                              rx=r"urine", name_search=r"(urine.*culture|culture.*urine)"),
}

# ---------------------------------------------------------------------------------------------------------
# DRUGS: RxNorm ingredient names (P§6.1, Appendix A.3), systemic routes only.
# ---------------------------------------------------------------------------------------------------------
ABX_CLASSES = {
    "A": ["nitrofurantoin", "fosfomycin"],
    "B": ["trimethoprim", "sulfamethoxazole"],
    "C": ["ciprofloxacin", "levofloxacin", "ofloxacin", "moxifloxacin", "norfloxacin", "gemifloxacin", "delafloxacin"],
    "D1": ["amoxicillin", "clavulanate", "ampicillin", "sulbactam", "penicillin V", "penicillin G", "dicloxacillin",
           "piperacillin", "tazobactam", "amdinocillin pivoxil"],
    "D2": ["cephalexin", "cefadroxil", "cefazolin", "cefuroxime", "cefaclor", "cefprozil", "cefpodoxime", "cefdinir",
           "cefixime", "cefditoren", "ceftriaxone", "cefotaxime", "ceftazidime", "cefepime", "cefoxitin", "cefotetan",
           "ceftaroline", "cefiderocol"],
    "D3": ["meropenem", "ertapenem", "imipenem", "aztreonam", "sulopenem", "sulopenem etzadroxil"],
    "E": ["doxycycline", "minocycline", "tetracycline"],
    "F": ["azithromycin", "clarithromycin", "erythromycin"],
    "G": ["clindamycin", "metronidazole", "tinidazole", "gentamicin", "tobramycin", "amikacin", "plazomicin",
          "vancomycin", "linezolid", "gepotidacin"],
}
# Companion ingredients that never define a class on their own (a clavulanate-only product does not exist).
COMPANION_ONLY = {"clavulanate", "sulbactam", "tazobactam"}
# P§5.3 antibiotics with urinary activity (confirm a UTI episode). Classes A, B, C, D1, D2, D3 + aminoglycosides + gepotidacin.
URINARY_ACTIVE = (ABX_CLASSES["A"] + ABX_CLASSES["B"] + ABX_CLASSES["C"] + ABX_CLASSES["D1"] + ABX_CLASSES["D2"]
                  + ABX_CLASSES["D3"] + ["gentamicin", "tobramycin", "amikacin", "plazomicin", "gepotidacin"])
# P§6.1 "not counted as antibiotic courses"
PROPHYLAXIS_ONLY = ["methenamine"]
NONABSORBED = ["fidaxomicin", "rifaximin", "neomycin"]                 # + oral vancomycin (dose-form rule below)
# P§6.5 prophylaxis agents
PROPHYLAXIS_AGENTS = ["nitrofurantoin", "trimethoprim", "cephalexin", "methenamine"]   # TMP-SMX counted under trimethoprim
# P§10 negative-control exposures
NEG_CONTROL_RX = ["fluconazole", "terbinafine", "acyclovir", "valacyclovir", "oseltamivir"]
# P Appendix A.4 kidney-cancer systemic drugs (outcome confirmation; prior-cancer exclusion)
RCC_DRUGS = ["sunitinib", "pazopanib", "axitinib", "cabozantinib", "lenvatinib", "sorafenib", "tivozanib", "everolimus",
             "temsirolimus", "belzutifan", "pembrolizumab", "nivolumab", "ipilimumab", "avelumab"]
# Covariate drug groups by ATC code (proposal: the protocol names the groups, not their ingredients). Name-checked in F00.
ATC_GROUPS = {
    "GLUCOSE_LOWERING": dict(atc=["A10"], rx=r"diabetes"),
    "ANTIHYPERTENSIVE": dict(atc=["C02", "C03", "C07", "C08", "C09"], rx=r"antihypertensive|diuretic|beta blocking|calcium channel|renin"),
    "SYSTEMIC_STEROID": dict(atc=["H02AB"], rx=r"glucocorticoid"),
    "IMMUNOSUPPRESSANT": dict(atc=["L04"], rx=r"immunosuppress"),
}
# Dose forms that are not systemic exposure (P§6.1). Matched against the drug concept name.
NON_SYSTEMIC_FORM_RX = (r"ophthalmic|otic|topical|vaginal|inhal|nasal|cream|ointment|lotion|\bgel\b|shampoo|"
                        r"irrigation|eye drop|ear drop|mouthwash|rectal|enema|transdermal|medicated pad|powder for nebuliz")
ORAL_FORM_RX = r"\boral\b|tablet|capsule|chewable|suspension|syrup"    # oral vancomycin -> non-absorbed

# ---------------------------------------------------------------------------------------------------------
# SURVEY (The Basics, Lifestyle). Codes are PPI concept_code strings, NOT concept IDs. Every item is
# [UNCONFIRMED] until F01 prints the v9 inventory; F02 stops on a code not found. Fill CROSSWALKS from F01.
# ---------------------------------------------------------------------------------------------------------
SURVEY_SEARCH_RX = (r"race|ethnic|hispanic|latin|puerto|birthplace|born|country|sex at birth|biological sex|education|"
                    r"highest grade|income|insurance|health insurance|smok|cigarette|zip")
SURVEY_ITEMS = {        # question concept_code -> variable. [UNCONFIRMED] codes, check against the F01 inventory.
    "Race_WhatRaceEthnicity": "race_eth",
    "BiologicalSexAtBirth_SexAtBirth": "sex_at_birth",
    "EducationLevel_HighestGrade": "education",
    "Income_AnnualIncome": "income",
    "Smoking_100CigsLifetime": "smoke_100",
    "StreetAddress_PIIZIP": "zip3",
    "Insurance_HealthInsurance": "insured",                     # confirmed in F01.6 (2026-11-28)
    "HealthInsurance_InsuranceTypeUpdate": "insurance_type",    # multi-select; confirmed in F01.6
}
# Insurance type (descriptive, P§9.4). Multi-select: the first group in INSURANCE_PRECEDENCE that a person ticked wins.
INSURANCE_GROUPS = {"InsuranceTypeUpdate_Medicaid": "Medicaid", "InsuranceTypeUpdate_Medicare": "Medicare",
                    "InsuranceTypeUpdate_Military": "Military or VA", "InsuranceTypeUpdate_VA": "Military or VA",
                    "InsuranceTypeUpdate_EmployerOrUnion": "Private", "InsuranceTypeUpdate_Purchased": "Private",
                    "InsuranceTypeUpdate_Indian": "Other", "InsuranceTypeUpdate_OtherHealthPlan": "Other",
                    "InsuranceTypeUpdate_None": "None"}
INSURANCE_PRECEDENCE = ["Medicaid", "Medicare", "Military or VA", "Private", "Other", "None"]
HISPANIC_SUBGROUP_QUESTION_RX = r"hispanic|latino|spanish"     # P§11.1: subgroup item located by name in F01
BIRTHPLACE_QUESTION_RX = r"born|birthplace|country of birth"
CROSSWALKS = {          # answer concept_code -> level; unmatched answers are printed by F04 (never silently dropped)
    "race_eth": {"WhatRaceEthnicity_White": "White", "WhatRaceEthnicity_Black": "Black",
                 "WhatRaceEthnicity_Hispanic": "Hispanic", "WhatRaceEthnicity_Asian": "Asian",
                 "WhatRaceEthnicity_AIAN": "AIAN", "WhatRaceEthnicity_NHPI": "NHPI", "WhatRaceEthnicity_MENA": "MENA",
                 "WhatRaceEthnicity_RaceEthnicityNoneOfThese": "Other", "PMI_Skip": "Not answered",
                 "PMI_PreferNotToAnswer": "Not answered", "PMI_DontKnow": "Not answered"},
    "sex_at_birth": {"SexAtBirth_Male": "Male", "SexAtBirth_Female": "Female", "SexAtBirth_Intersex": "Other",
                     "SexAtBirth_SexAtBirthNoneOfThese": "Other", "PMI_Skip": "Not answered",
                     "PMI_PreferNotToAnswer": "Not answered", "PMI_DontKnow": "Not answered"},
    # Answer codes confirmed in the 2026-11-28 F01.6 inventory; groupings approved 2026-11-28. Every survey question:
    # skipped, prefer not to answer and don't know = "Not answered".
    "education": {"HighestGrade_NeverAttended": "Less than high school", "HighestGrade_OneThroughFour": "Less than high school",
                  "HighestGrade_FiveThroughEight": "Less than high school",
                  "HighestGrade_NineThroughEleven": "Less than high school",
                  "HighestGrade_TwelveOrGED": "High school or GED", "HighestGrade_CollegeOnetoThree": "Some college",
                  "HighestGrade_CollegeGraduate": "College graduate", "HighestGrade_AdvancedDegree": "Advanced degree",
                  "PMI_Skip": "Not answered", "PMI_PreferNotToAnswer": "Not answered", "PMI_DontKnow": "Not answered"},
    "income": {"AnnualIncome_less10k": "<25k", "AnnualIncome_10k25k": "<25k", "AnnualIncome_25k35k": "25-50k",
               "AnnualIncome_35k50k": "25-50k", "AnnualIncome_50k75k": "50-100k", "AnnualIncome_75k100k": "50-100k",
               "AnnualIncome_100k150k": "100-200k", "AnnualIncome_150k200k": "100-200k", "AnnualIncome_more200k": ">=200k",
               "PMI_Skip": "Not answered", "PMI_PreferNotToAnswer": "Not answered", "PMI_DontKnow": "Not answered"},
}
SKIP_CODES = ["PMI_Skip", "PMI_PreferNotToAnswer", "PMI_DontKnow"]     # kept as "Not answered" (M9)
PR_ZIP3 = ["006", "007", "009"]                                        # P§11.1 [UNCONFIRMED against USPS]

# ---------------------------------------------------------------------------------------------------------
# ADJUSTMENT SETS (P§9.2). Variable names are the columns F04 builds.
# ---------------------------------------------------------------------------------------------------------
ADJ_Q1 = ["age_t0", "sex", "t0_year", "smoking", "bmi_cat", "hypertension", "diabetes", "ckd_pre_e1", "stones",
          "immunosuppression", "visits_prior_year", "imaging_prior_3y", "prior_uti_episodes", "prior_abx_courses",
          "education", "income", "race_eth"]                 # + site as stratum (P§9.2)
ADJ_Q2 = ADJ_Q1 + ["episode_type", "episode_setting", "resistant_org", "pcn_allergy", "egfr_cat"]

TABLE1_VARS = ADJ_Q2 + ["coding_era", "years_t0_to_enrollment", "ehr_years_before_t0", "hispanic", "puerto_rican",
                        "pr_residence", "insurance", "prophylaxis_before", "hematuria_2y"]
# Table 1 standardized differences (blind: covariates only): (group, reference) pairs
SMD_PAIRS = {"Q1": [("recurrent", "single")], "Q2": [("B", "A"), ("C", "A"), ("D", "A")]}

# ---------------------------------------------------------------------------------------------------------
# RUN SETTINGS
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# OUTCOME MODELS (F05-F10). Choices the protocol leaves to the analyst are README open items O15-O19.
# ---------------------------------------------------------------------------------------------------------
ANALYSIS = {
    "min_events_per_arm": DECISIONS["D15_subgroup_min_per_cell"],   # an HR is shown only if each compared arm has >= 21
    "ps_truncate_pct": (1, 99),       # O15: stabilized IPTW truncated at the 1st and 99th percentiles
    "site_stratum": True,             # P§9.2: EHR site as Cox stratum in the primary model (A14 shows it unstratified)
    "mice_m": 20,                     # P§12.2 #8 / D16: primary
    "mice_m_sensitivity": 5,          # O16: imputations for sensitivity and secondary models (compute)
    "mice_iters": 5,
    "mice_fit_rows": 50000,           # O16: rows sampled to fit each imputation model (donors are all observed rows)
    "cif_times": (5.0, 10.0),         # years after outcome counting starts (P§2.1)
    "cif_bootstrap": 100,             # per imputation
    "ph_split_years": 3.0,            # period HRs 2-5 and > 5 years after T0 = 0-3 and > 3 years after S (P§12.2 #6)
    "count_covariates_log1p": ["visits_prior_year", "imaging_prior_3y", "prior_uti_episodes", "prior_abx_courses"],
    "q2b_course_bins": [0, 1, 2, 6],  # D10: 0 / 1 / 2-5 / >= 6 cumulative courses
    "well_observed_history_years": 5,
}
ORDINAL_LEVELS = {     # imputation order (P§9.4 groupings approved 2026-11-28)
    "education": ["Less than high school", "High school or GED", "Some college", "College graduate", "Advanced degree"],
    "income": ["<25k", "25-50k", "50-100k", "100-200k", ">=200k"],
}

RUN = {
    # Freeze gate (STUDY-LIFECYCLE stage 8; P§12.1): file name in workbench_outputs/ of the blind F00-F04 output.
    # F05-F10 stop until it is set, the file is in the repo and its F04 code hash matches this run's F04.
    "blind_frozen_output": "",
    # PI deviation 2026-11-28 (docs/AMENDMENTS-RUTI-RCC-2026-11-28.md, DEV1): the blind freeze is waived and F05-F10
    # run on the same run folder's F04 cohorts. Non-empty = waived; the text is printed in every outcome cell.
    "blinding_waived": "DEV1 2026-11-28: PI waived the blind freeze after the D17 gate (Q1 35, Q2 30 events); "
                       "all outcome estimates are exploratory (D17 < 100 events)",
    "max_gib": 150,                # EHR span scans ~76 GiB in v9 (2026-11-28 run); about USD 0.50 per scan on demand                 # refuse a query that would scan more (q(max_gib=))
    "bucket_name": "",             # workspace bucket resource id; run folder goes there (WORKBENCH-2.0-GUIDE §4)
    "cdr_resource": "",            # data collection resource id for `wb resolve` (WORKBENCH-2.0-GUIDE §2)
}
