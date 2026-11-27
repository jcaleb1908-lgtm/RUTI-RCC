##### BEGIN F00 setup and concept sets #####
# Preflight, config snapshot, every code set resolved by code / ingredient name inside the CDR, the standard
# concepts reached in each coding era, and the Appendix D check sheet. Saves concept tables (not participant-level).
import json, hashlib
import pandas as pd
import aou_core as C
import _config as K
import rutirc as R

C.CFG.update(project=K.PROJECT, version=K.PIPELINE_VERSION, hide_max=20)
C.begin("F00", "setup and concept sets", create=RUN_CONTROL.get("new_run", True))
OK = C.bq_preflight()
if not OK:
    raise SystemExit("STOPPED: BigQuery preflight failed (see message above).")
if K.CDR_EXPECTED not in C.CFG["cdr"]:
    print(f"WARNING: CFG['cdr'] = {C.CFG['cdr']} does not contain {K.CDR_EXPECTED} (P§14). Confirm the workspace CDR.")


def p_snapshot():
    """Config snapshot and code hash: every later output can be traced to these definitions."""
    snap = {k: getattr(K, k) for k in ["PIPELINE_VERSION", "PROTOCOL_VERSION", "DECISIONS", "IMPLEMENTATION", "AGE_MIN"]}
    txt = json.dumps(snap, default=str, indent=1)
    (C._run() / "config_snapshot.json").write_text(txt)
    print("code hash:", CODE_HASH)
    print(txt[:3000])
C.part("1 config snapshot", p_snapshot)


def p_codes():
    """Resolve every ICD-9-CM / ICD-10-CM / CPT4 / ICD-10-PCS / ICD-9 procedure / LOINC / ATC code by code."""
    con = C.q(R.CONCEPT_BY_CODE_SQL, [C.A("vocabs", "STRING", R.CODE_VOCABS), C.S("rx", "STRING", R.all_code_regex())],
              "concepts by code", max_gib=K.RUN["max_gib"])
    S = R.assign_code_sets(con)
    C.save_state(S, "code_sets", internal=False)
    cov = R.code_set_coverage(S)
    C.show("1", "Protocol codes not found in this CDR (codes_found = 0)", cov.loc[cov["codes_found"].eq(0)],
           counts=[], note="A code absent here is never matched; check it is not a typo before the run.")
    bad = S.loc[~S["name_ok"], ["set_key", "vocabulary_id", "concept_code", "concept_name"]]
    C.show("2", "Codes whose name does not match the set's expected-name regex (read each one)", bad, counts=[])
    C.show("3", "Codes per set and vocabulary", S.groupby(["set_key", "vocabulary_id"]).size().rename("codes").reset_index(),
           counts=[])
C.part("2 code sets", p_codes)


def p_maps():
    """Standard concepts the source codes map to, by era (Appendix D 'standard concepts reached'; M3)."""
    S = C.load_state("code_sets", internal=False)
    src = S.loc[S["vocabulary_id"].isin(["ICD9CM", "ICD10CM", "ICD9Proc", "ICD10PCS"])]
    M = C.q(R.MAPS_TO_SQL, [C.A("ids", "INT64", src["concept_id"].astype(int).unique().tolist())], "maps to")
    M = M.merge(src[["concept_id", "set_key", "era", "concept_code"]], left_on="source_concept_id", right_on="concept_id")
    C.save_state(M, "maps_to", internal=False)
    reach = M.groupby(["set_key", "era", "standard_name", "standard_domain"]).size().rename("source_codes").reset_index()
    C.show("4", "Standard concepts reached per set and era (read by name)", reach, counts=[], max_rows=400)
    unm = src.loc[~src["concept_id"].isin(M["source_concept_id"])]
    C.show("5", "Source codes with no 'Maps to' (still matched by source code)", unm[["set_key", "era", "concept_code", "concept_name"]],
           counts=[])
C.part("3 source-to-standard map", p_maps)


def p_drugs():
    """RxNorm ingredients by name, all descendants, dose-form class and study groups (M5: read routes)."""
    names = sorted({n.lower() for v in K.ABX_CLASSES.values() for n in v} | set(K.PROPHYLAXIS_ONLY) | set(K.NONABSORBED)
                   | set(K.NEG_CONTROL_RX) | set(K.RCC_DRUGS))
    ing = C.q(R.INGREDIENT_SQL, [C.A("names", "STRING", names)], "ingredients by name")
    missing = sorted(set(names) - set(ing["concept_name"].str.lower()))
    C.show("6", "RxNorm ingredient names not found (P Appendix A.3)", pd.DataFrame({"ingredient": missing}), counts=[])
    desc = C.q(R.DRUG_DESC_SQL, [C.A("ids", "INT64", ing["concept_id"].astype(int).tolist())], "drug descendants",
               max_gib=K.RUN["max_gib"])
    G = R.drug_glossary(ing, desc)
    C.save_state(G, "drug_glossary", internal=False)
    t = G.groupby(["abx_classes", "form", "abx_counted"]).size().rename("drug_concepts").reset_index()
    C.show("7", "Antibiotic drug concepts by class, dose form and whether counted", t, counts=[])
    C.show("8", "Products in two antibiotic classes (read; counted in each)", G.loc[G["multi_class"], ["concept_name", "abx_classes"]],
           counts=[])
    C.show("9", "Forms classified unknown (read; counted as systemic)",
           G.loc[G["is_abx"] & G["form"].eq("unknown_form"), ["concept_name", "abx_classes"]].head(200), counts=[])
    # ATC covariate groups -> drug descendants
    S = C.load_state("code_sets", internal=False)
    atc = S.loc[S["kind"].eq("ATC")]
    A = C.q(R.DRUG_DESC_SQL, [C.A("ids", "INT64", atc["concept_id"].astype(int).tolist())], "ATC descendants")
    A = A.merge(atc[["concept_id", "set_key"]].rename(columns={"concept_id": "root_id"}), on="root_id")
    C.save_state(A[["set_key", "concept_id", "concept_name", "concept_class_id"]].drop_duplicates(), "atc_drugs", internal=False)
    C.show("10", "ATC covariate groups: drug concepts reached", A.groupby("set_key").size().rename("drug_concepts").reset_index(),
           counts=[])
C.part("4 drug concept sets", p_drugs)


def p_appendix_d():
    """Appendix D sheet, one row per definition; builder and checker initials are filled by hand."""
    S = C.load_state("code_sets", internal=False)
    rows = []
    for key, g in S.groupby("set_key"):
        rows.append({"definition": f"{key} {K.PIPELINE_VERSION}", "starting_codes": len(g),
                     "icd9_codes": int(g["era"].eq("ICD-9").sum()), "icd10_codes": int(g["era"].eq("ICD-10").sum()),
                     "name_check_failures": int((~g["name_ok"]).sum()), "wrong_concept_removed": "", "checked_by": ""})
    C.show("11", "Appendix D concept-set check sheet (fill the last two columns by hand)", pd.DataFrame(rows), counts=[])
C.part("5 Appendix D sheet", p_appendix_d)
C.end()
##### END F00 #####
