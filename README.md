# RUTI-RCC: recurrent UTI, its antibiotics, and kidney cancer in All of Us

**Question.** Is the kidney cancer (RCC) risk seen after recurrent UTI explained by the infection, by the antibiotics used to treat it, by both, or by detection and reverse causation? Does it differ in Hispanic and Puerto Rican participants?
**Data.** All of Us Controlled Tier, CDR v9 **C2025Q4R6** (data through 2025-01-01), Researcher Workbench 2.0.
**Protocol.** [`docs/PROTOCOL-RUTI-RCC-2026-11-19-v0.1.md`](docs/PROTOCOL-RUTI-RCC-2026-11-19-v0.1.md), a draft for PI review that is **not signed yet**. The literature review v3 is in `docs/`.
**Team.** PI: Dra. Yaritza Inostroza-Nieves. Protocol: Valeria Román. Workspace and code: Jonathan.
**Lifecycle stage** (BRAIN `projects/all-of-us/STUDY-LIFECYCLE.md`): 4 (protocol) → 6 (code). The code reaches the blind line. Stages 7 and 8 (the cohort description run and the feasibility gate) can run now. Outcome models wait for the PI signature.

This README is the only notes file: facts, decisions (D#), open items (O#) and limitations (L#).

---

## 1. What is built (v0.1, 2026-11-27)

| Cell | Does | Protocol |
|---|---|---|
| `F00_setup` | Preflight and config snapshot. Resolves every code set **by code or ingredient name** inside the CDR (no typed concept IDs). Prints names that fail the name check, protocol codes missing from the CDR, the standard concepts each ICD-9 and ICD-10 code maps to, drug concepts by class and dose form, and the Appendix D check sheet | App. A, C, D; §6.1 |
| `F01_inventory` | Whole-CDR data checks (F7, §16.2): date-of-birth precision, src_id, days-supply and route completeness by drug type, urine culture by name, death sources, the Basics and Lifestyle survey inventory (including the Hispanic subgroup and Puerto Rico answers), PR ZIP3, site-year drug completeness, observation_period | §13 F7, §16.2 |
| `F02_extract` | Candidates = anyone with a §5.1 UTI code or an RCC code. Pulls diagnoses from 4 tables, procedures, drugs, visits, EHR span, activity months, labs, person, death, enrollment, survey | §3–9 |
| `F03_phenotypes` | UTI episodes (31-day fixed window; code + antibiotic ±3 d or + culture ±7 d), antibiotic courses (7-day merge, one per inpatient stay, ±3 d indication), prophylaxis (3 rules, 90-day blocks), negative-control courses, eGFR (CKD-EPI 2021), ESRD and CKD dates, RCC confirmed / possible / ambiguous with diagnosis and detection dates | §4, 5, 6, 8 |
| `F04_cohorts` | One `build(spec)` for Q1 (landmark at E1 + 365 d) and Q2 (new-user active comparator). Masked flow, Table 1 by exposure, arm sizes, timing. **RCC events in total only (blind)**. D17 go/no-go and events-per-parameter gate | §2, 3, 7, 12.1, 13 |

The cohort builder follows the protocol's order: candidates, T0 window, age ≥ 40, 3-year history with a visit in each year, data availability, exclusions, no kidney cancer before T0, observed at T0, Q2 washout, one index per person, then the 2-year latency, then follow-up (RCC, death, last EHR activity, ESRD and the CDR cutoff) and covariates measured before E1 (Q1) or T0 (Q2). The index is chosen **before** the latency step, so a person's index never depends on a later outcome.

**Not built yet (phase 2, after the PI signs):** F05 IPTW Cox plus Aalen–Johansen and PH checks; F06 lag/latency variants, delayed entry, Q2b time-varying and Q2c; F07 MICE (m = 20); F08 sex and Hispanic/PR modification with RERI; F09 sensitivity, negative controls, E-values, site and validation; F10–F11 glossary and tables. `docs/ANALYSIS-MAP-RUTI-RCC-2026-11-27.csv` maps all 22 protocol analyses to their cells.

## 2. How to run (Workbench 2.0)

1. **Locally, first:** `pip install -r requirements.txt && python pipeline/test_synthetic.py` should print `all tests passed`. The test builds a synthetic OMOP CDR in DuckDB, runs the real BigQuery SQL through sqlglot, runs F00–F04 end to end, and scans every saved table for an unmasked count of 1–20.
2. In the workspace: add this repo as a Git resource (Apps → Git repositories). Add the v9 data collection. Create a bucket for runs.
3. In a JupyterLab app terminal: `wb resource list`, then `wb resolve --id=<cdr id>`. Raise autostop before the run.
4. Paste `pipeline/RUN_IN_WORKBENCH.py` into a notebook cell. Set `REPO`, `CDR` and `BUCKET` at the top, then run it. After a kernel restart, `RUN_CONTROL["resume_from"] = "F03"` resumes from that cell.
5. **Run F00 and F01 alone the first time** (`run_only=["F00_setup", "F01_inventory"]`). Read every name-check table and the survey inventory, then fill §4 below and `_config.py` (`SURVEY_ITEMS`, `CROSSWALKS`) before running F02–F04.
6. Copy out **only** `runs/<run>/tables/ALL_TABLES.txt` and save it as `workbench_outputs/YYYY-MM-DD_<what>_output.txt`. Add a row to the run log (§6).

Data rules: participant-level data never leave the Workbench, and no AI tool runs inside the workspace (DUCC; WORKBENCH-2.0-GUIDE §7). Every table goes through `aou_core.show()` (counts 1–20 → `<=20`, complementary cells hidden). The flow goes through `flow_mask`. A second person checks each output before it is shared (P§14).

## 3. Decisions (PI, protocol §16)

Every D1–D21 decision is a named switch in `pipeline/_config.py` `DECISIONS`, set to the protocol's **proposed** value. None is PI-approved yet. **D18 is urgent:** Jonathan's free credits expire **2026-12-24**, and buckets may be deleted when they do. Don't build cohorts on that workspace; move to a PI- or grant-funded pod first.

| D# | Default in code | Status |
|---|---|---|
| D1 one paper | Q1, Q1b, Q2, Q2b main; Q3 exploratory | PI |
| D2–D5 | landmark 365 d; ref. nitrofurantoin/fosfomycin; T0 before enrollment allowed; visit in each of 3 years; T0 2008–2022; 365-d clean window | PI |
| D6–D10 | hereditary excluded only if coded before T0; 31-d episode; code + abx ±3 d or culture ±7 d; N39.0 = lower; any-type recurrence primary; courses 0/1/2–5/≥6 | PI |
| D11–D17, D20 | 2-y lag and latency; confirmed RCC, C65 excluded; both sexes; NCOs adrenal adenoma and cataract; subgroup ≥100 events; MICE; 300/100 gate; ESRD censors | PI |
| D18, D19, D21 | funding, IRB, journal and OSF | PI (not code) |

## 4. Open items (O#): choices the protocol leaves open, implemented with a default

| O# | Item | Default now | Who settles |
|---|---|---|---|
| O1 | "First qualifying UTI episode": the person's first ever, or the first whose T0 passes eligibility? | first **eligible** (`O1_index_choice`) | PI |
| O2 | Does an *unconfirmed* coded episode in the prior 365 d block E1? | yes (conservative) | PI |
| O3 | "≥2 episodes within 180 days": ≤ 180 or < 180 | ≤ 180 | PI |
| O4 | Q2 "≥2 classes within ±3 d": among A–G or only the compared classes | A–G | PI |
| O5 | Which src_id values count as EHR (`^EHR site`) | **Checked 2026-11-28 (F01.2):** `EHR site #` = 69 sites, counted. "Participant Mediated EHR" (PTSC, TPC, CLAD) and the portals are **not** counted as EHR. | PI to confirm that participant-mediated EHR stays out |
| O6 | Date of birth | **Settled 2026-11-28 (F01.1):** only the birth **year** is real. `day_of_birth` is empty for everyone, and `birth_datetime` has one month and day 15 for all 747,029 persons (a fixed placeholder). Age is therefore accurate to about ±1 year (L11) | none: data limit |
| O7 | RCC diagnosis date = earliest RCC code (literal §8.4), not the code that met confirmation | literal | PI |
| O8 | ≥2 RCC codes that don't meet 30–365 d (e.g. 4 d apart) | "possible" | PI |
| O9 | Codes the protocol named by group only: antihypertensive / glucose-lowering / steroid / immunosuppressant by ATC (C02–C09, A10, H02AB, L04); resistant organism Z16 / V09; transplant status Z94.0, T86.1, V42.0, 996.81 (invalidates confirming procedures); LOINC creatinine 2160-0 and 38483-4, HbA1c 4548-4 and 17856-6, BMI 39156-5 | **Approved 2026-11-28**; written into protocol A.1 and A.5 | done |
| O10 | Survey PPI codes and groupings | **Settled 2026-11-28:** all `SURVEY_ITEMS` codes exist; Puerto Rican subgroup released. Education (5 levels) and income (5 brackets) groupings approved. Every question: skip / prefer not / don't know = "Not answered" (protocol 9.4) | done |
| O11 | Person-level "data available" rule measured by month (≥1 EHR condition and ≥1 EHR drug record in the 36 months before T0) | month precision | note |
| O12 | People whose ESRD date falls before outcome counting starts are left out of the analysed cohort, like death | yes | PI |
| O13 | Codes F00.2 flagged | **Settled 2026-11-28:** keep E66.3 (moot) and N02.B (IgA nephropathy causes hematuria); **drop 366.11** (cataract NCO = only codes whose ICD title is a cataract); **keep nephroureterectomy** (CPT 50234, 50236, 50548; ICD-9 55.51) as confirming, with a sensitivity without it: F04 reports the events that depend on it (protocol 8.3) | done |
| O14 | Urine culture by LOINC **and** SNOMED (codes in protocol 5.3; code systems explained in Appendix A) | **Approved 2026-11-28** | done |

## 5. Limitations (candidates, L#)

L1 volunteer cohort (AL1) · L2 partial EHR capture: outpatient antibiotics and UTIs are under-counted (AL2; P§3.3) · L3 survey covariates collected at enrollment, after T0 (AL4) · L4 left truncation when T0 is before enrollment (AL5; D4) · L5 site heterogeneity (AL7) · L6 death under-ascertained (AL8) · L7 no validated US EHR RCC algorithm; codes carry no histology (AL9; P§8.1) · L8 ICD-9 mapping gaps (AL10; F00 table 4) · L9 hereditary RCC syndromes without specific codes cannot be excluded (P§4) · L10 nitrofurantoin is not a pure no-drug reference (NTP TR-341; P§1.1) · L11 birth year only (O6): age at T0 can be off by up to a year, which blurs the age ≥ 40 cut-off · L12 drug capture varies by site: 309 of 1,554 EHR site-years have visits but no drug records (F01.11).

## 6. Run log

| Run ID | Date (UTC) | Code hash | Commit | CDR | Cells | Ended | Output | Kept? |
|---|---|---|---|---|---|---|---|---|
| first F00+F01 | 2026-11-28 | `151a1d6eb424` | `872100a` | C2025Q4R6 | F00, F01 | ok | `workbench_outputs/2026-11-28_F00-F01_output.txt` | yes |
| F00+F01 rerun | 2026-11-28 | `a58e5386e92a` | `d368550` | C2025Q4R6 | F00, F01 | ok, but `_config`/`rutirc` were still the first run's (kernel not restarted) | `workbench_outputs/2026-11-28_F00-F01_rerun_output.txt` | partial: F00.12 and F01.11 only |
| F00+F01 final | 2026-11-28 | `2485c8c8ed67` | `c74c4dc` (PR #4) | C2025Q4R6 | F00, F01 | ok | `workbench_outputs/2026-11-28_F00-F01_final_output.txt` | yes: data checks complete |

## 7. Provenance

- `pipeline/aou_core.py` is copied unchanged from BRAIN `projects/all-of-us/starter/aou_core.py` at BRAIN commit `035e4d7` (sha256 `492cdc51…`), plus the 2026-11-28 `bucket_mount` fix (no recursive walk of `~/workspace`; sha256 `7609b5c1…`), applied to BRAIN too. A fix there is a fix here too (M39).
- Amendments after the PI signs go in `docs/AMENDMENTS-RUTI-RCC-<date>.md` (BRAIN template).
