# RUTI-RCC: recurrent UTI, its antibiotics, and kidney cancer in All of Us

**Question.** Is the kidney cancer (RCC) risk seen after recurrent UTI explained by the infection, by the antibiotics used to treat it, by both, or by detection and reverse causation? Does it differ in Hispanic and Puerto Rican participants?
**Data.** All of Us Controlled Tier, CDR v9 **C2025Q4R6** (data through 2025-01-01), Researcher Workbench 2.0.
**Protocol.** [`docs/PROTOCOL-RUTI-RCC-2026-11-19-v0.1.md`](docs/PROTOCOL-RUTI-RCC-2026-11-19-v0.1.md), **D1–D21 approved by the PI on 2026-11-28, all as proposed** (reported by Jonathan). The PI's written signature on the document is still to be added; the commit that records this approval is the prespecification stamp. The literature review v3 is in `docs/`.
**Team.** PI: Dra. Yaritza Inostroza-Nieves. Protocol: Valeria Román. Workspace and code: Jonathan.
**Lifecycle stage** (BRAIN `projects/all-of-us/STUDY-LIFECYCLE.md`): 6 (code), protocol decisions approved. Next: stages 7 and 8 (the blind cohort description run, F02–F04, and the feasibility gate). The outcome models (F05–F10) are written and refuse to run until that blind output is saved in `workbench_outputs/` and named in `_config.RUN["blind_frozen_output"]`.

This README is the only notes file: facts, decisions (D#), open items (O#) and limitations (L#).

---

## 1. What is built (v0.2, 2026-11-28)

| Cell | Does | Protocol |
|---|---|---|
| `F00_setup` | Preflight and config snapshot. Resolves every code set **by code or ingredient name** inside the CDR (no typed concept IDs). Prints names that fail the name check, protocol codes missing from the CDR, the standard concepts each ICD-9 and ICD-10 code maps to, drug concepts by class and dose form, and the Appendix D check sheet | App. A, C, D; §6.1 |
| `F01_inventory` | Whole-CDR data checks (F7, §16.2): date-of-birth precision, src_id, days-supply and route completeness by drug type, urine culture by name, death sources, the Basics and Lifestyle survey inventory (including the Hispanic subgroup and Puerto Rico answers), PR ZIP3, site-year drug completeness, observation_period | §13 F7, §16.2 |
| `F02_extract` | Candidates = anyone with a §5.1 UTI code or an RCC code. Pulls diagnoses from 4 tables, procedures, drugs, visits, EHR span, activity months, labs, person, death, enrollment, survey | §3–9 |
| `F03_phenotypes` | UTI episodes (31-day fixed window; code + antibiotic ±3 d or + culture ±7 d), antibiotic courses (7-day merge, one per inpatient stay, ±3 d indication), prophylaxis (3 rules, 90-day blocks), negative-control courses, eGFR (CKD-EPI 2021), ESRD and CKD dates, RCC confirmed / possible / ambiguous with diagnosis and detection dates | §4, 5, 6, 8 |
| `F04_cohorts` | One `build(spec)` for Q1 (landmark at E1 + 365 d) and Q2 (new-user active comparator). Masked flow, Table 1 by exposure, arm sizes, timing. **RCC events in total only (blind)**. D17 go/no-go and events-per-parameter gate | §2, 3, 7, 12.1, 13 |

| `F05_main` | Rates per 1,000 PY by group (A3); IPTW Aalen–Johansen cumulative incidence at 5 and 10 y and the difference (A4); **primary IPTW cause-specific Cox, MI m = 20, site stratum, robust SE** (A5, A5b); PH check: score test on log t and period HRs 2–5 / > 5 y after T0 (A6); D9 UTI types, Q1b, Q2 head-to-head (A7a–c) | §2, 12.2 #3–7 |
| `F06_timevarying` | Q2b: time-varying cumulative courses by class (0/1/2–5/≥6, lagged 2 y) plus a per-course trend, adjusted for ADJ_Q1 and lagged UTI episodes (A7d); latency 1 and 3 y and the Q2b lag × latency grid (A9); delayed entry at enrollment (A13); Q2c marked not run (A21) | §2.3, 12.2 #7, 9, 13 |
| `F07_missing` | Complete case and missing-as-category versions of the primary contrasts; imputation diagnostics (A8) | §12.2 #8, D16 |
| `F08_modification` | Sex-specific HRs and exposure × sex interaction with RERI (A19); Hispanic and Puerto Rican group HRs, product-term and RERI interaction and 10-y CIF difference where the D15 gate passes, rates by group otherwise, with a cross-table leak check (A20) | §9.3, 11, 12.2 #19–20 |
| `F09_sensitivity` | Confirmed + possible, detection-dated, without nephroureterectomy (A10); 14-day window, code-only, culture-only UTI (A11); well observed (A12); site unstratified, leave-one-site-out, site-year completeness (A14); ICD-10 era (A15); no stones, no hematuria, outpatient index antibiotic, censor at new cancer, Fine–Gray (A16); time-varying post-T0 imaging and by UTI type (A17); negative-control outcomes and exposure (A18); validation of confirmation routes (A22) | §10, 12.2 #10–18, 22 |
| `F10_results` | Every HR in one table with E-values; the §12.1 robustness rule for each primary contrast; glossary | §10, 12.1 |

Estimation code: `pipeline/analysis.py` (weighted Cox with entry, strata and cluster-robust SE; IPTW; balance; Aalen–Johansen; counting-process splits; Fine–Gray weighting; MICE with predictive mean matching; RERI; Rubin's rules) and `pipeline/outcomes.py` (the one workflow every outcome cell calls). The Cox code is tested against statsmodels PHReg (unweighted) and against row replication (weights).

**Freeze gate.** F05–F10 stop unless the blind F00–F04 output is saved in `workbench_outputs/`, its file name is set in `_config.RUN["blind_frozen_output"]`, and its F04 code hash is the F04 that built the current run's cohorts (`state/F04_stamp.json`). So the models can only run on cohorts whose blind description was reviewed and committed first (STUDY-LIFECYCLE stage 8; P§12.1).

**Masking in the outcome cells.** An HR, CIF or RERI is shown only when each compared arm (or cell) has ≥ 21 events (O19). Event counts by arm appear once, in `F05.rates`; the sensitivity tables carry no counts, so nested analyses cannot be subtracted from each other. The F08 rates by ethnic group use `aou_core.leak_hidden` against the F05 arm totals.

**Not run, and why:** A21 Q2c needs a cohort of people with no UTI ever, which F02 does not extract. The survey part of A22 needs the Personal Medical History items, also not extracted. The ESRD inverse-probability-of-censoring sensitivity (protocol §4 table) is not written. Each can be added on request.

The cohort builder follows the protocol's order: candidates, T0 window, age ≥ 40, 3-year history with a visit in each year, data availability, exclusions, no kidney cancer before T0, observed at T0, Q2 washout, one index per person, then the 2-year latency, then follow-up (RCC, death, last EHR activity, ESRD and the CDR cutoff) and covariates measured before E1 (Q1) or T0 (Q2). The index is chosen **before** the latency step, so a person's index never depends on a later outcome.

`docs/ANALYSIS-MAP-RUTI-RCC-2026-11-27.csv` maps all 22 protocol analyses to their cells.

## 2. How to run (Workbench 2.0)

1. **Locally, first:** `pip install -r requirements.txt && python pipeline/test_synthetic.py` should print `all tests passed`. The test builds a synthetic OMOP CDR in DuckDB, runs the real BigQuery SQL through sqlglot, runs F00–F04 end to end, checks the freeze gate, runs F05–F10 end to end (with the ≥ 21-events rule lowered only for the tiny synthetic cohort), and scans every saved table for an unmasked count of 1–20.
2. In the workspace: add this repo as a Git resource (Apps → Git repositories). Add the v9 data collection. Create a bucket for runs.
3. In a JupyterLab app terminal: `wb resource list`, then `wb resolve --id=<cdr id>`. Raise autostop before the run.
4. Paste `pipeline/RUN_IN_WORKBENCH.py` into a notebook cell. Set `REPO`, `CDR` and `BUCKET` at the top, then run it. After a kernel restart, `RUN_CONTROL["resume_from"] = "F03"` resumes from that cell.
5. **Run F00 and F01 alone the first time** (`run_only=["F00_setup", "F01_inventory"]`). Read every name-check table and the survey inventory, then fill §4 below and `_config.py` (`SURVEY_ITEMS`, `CROSSWALKS`) before running F02–F04.
6. **Blind run:** `STAGE = "blind"` runs F00–F04 in a new run folder. Copy out ALL_TABLES.txt; it is reviewed and committed as the frozen blind output, and its name goes into `_config.RUN["blind_frozen_output"]`.
7. **Model run:** after `git pull` in the repo folder, set `STAGE = "models"`. F05–F10 run in the **same** run folder (the latest one). F09 is the longest cell; after a kernel restart use `RUN_CONTROL["resume_from"] = "F07"` and so on.
8. Copy out **only** `runs/<run>/tables/ALL_TABLES.txt` and save it as `workbench_outputs/YYYY-MM-DD_<what>_output.txt`. Add a row to the run log (§6).

Data rules: participant-level data never leave the Workbench, and no AI tool runs inside the workspace (DUCC; WORKBENCH-2.0-GUIDE §7). Every table goes through `aou_core.show()` (counts 1–20 → `<=20`, complementary cells hidden). The flow goes through `flow_mask`. A second person checks each output before it is shared (P§14).

## 3. Decisions (PI, protocol §16)

Every D1–D21 decision is a named switch in `pipeline/_config.py` `DECISIONS`, set to the protocol's proposed value. **The PI approved all of them as proposed on 2026-11-28.** D18 (funding): the PI expects funding if the analysis is significant, so the team runs the blind cohort on Jonathan's free credits, which expire **2026-12-24** (buckets may be deleted then). Copy each ALL_TABLES.txt out the day it runs, and move to a funded pod before the credits end.

| D# | Default in code | Status |
|---|---|---|
| D1 one paper | Q1, Q1b, Q2, Q2b main; Q3 exploratory | PI approved 2026-11-28 |
| D2–D5 | landmark 365 d; ref. nitrofurantoin/fosfomycin; T0 before enrollment allowed; visit in each of 3 years; T0 2008–2022; 365-d clean window | PI approved 2026-11-28 |
| D6–D10 | hereditary excluded only if coded before T0; 31-d episode; code + abx ±3 d or culture ±7 d; N39.0 = lower; any-type recurrence primary; courses 0/1/2–5/≥6 | PI approved 2026-11-28 |
| D11–D17, D20 | 2-y lag and latency; confirmed RCC, C65 excluded; both sexes; NCOs adrenal adenoma and cataract; subgroup ≥100 events; MICE; 300/100 gate; ESRD censors | PI approved 2026-11-28 |
| D18, D19, D21 | funding, IRB, journal and OSF | PI approved 2026-11-28 (not code); D18 see above |

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

| O15 | Propensity model and weights | Logistic (Q1) or multinomial logistic (Q2) on the §9.2 set, linear terms, log(1 + x) for the 4 count covariates; stabilized ATE weights truncated at the 1st and 99th percentiles | PI |
| O16 | Multiple imputation details | Predictive mean matching (5 donors, 5 cycles) for BMI, eGFR (Q2), education and income (ordinal); race "not answered" stays a group (§11.2) and is not imputed; m = 20 for A5/A5b and the CIF, m = 5 for secondary and sensitivity models (compute); each imputation model fitted on up to 50,000 sampled rows | PI |
| O17 | Regression-adjusted models (Q2b) with a covariate level that has no events, or exactly collinear columns | the level is merged into the most frequent level (a yes/no covariate leaves the model); collinear columns are dropped; top course categories with no events merge downward | note |
| O18 | EHR site | src_id of the last EHR visit before T0. §9.2 makes site a stratum (primary) while §12.2 #14 lists "stratified by site" as a sensitivity: the code stratifies the primary and shows the unstratified model in A14 | PI |
| O19 | When an estimate is shown | ≥ 21 events in each compared arm for every HR (the D15 per-cell rule applied to all analyses); event counts by arm only in F05.rates | note |
| O20 | Small definitions inside the analyses | Q2b counts courses and episodes from E1; CIF times are years after outcome counting starts; "well observed" = ≥ 5 y history and a visit in every complete follow-up year; the hematuria exclusion adds UTI episodes with a hematuria code; the PH test is an approximate score test on log t | PI |
| O21 | D17 result of the first complete blind run (2026-11-28): Q1 35 events (recurrent + single), Q2 30 events; events per parameter 0.9 and 0.6. Under D17 as approved, both are "descriptive + Q2b only" | F05–F10 code runs, but the primary HRs are not confirmatory under D17 | PI: follow D17, or amend the design before any outcome-by-exposure output (amendment file) |

## 5. Limitations (candidates, L#)

L1 volunteer cohort (AL1) · L2 partial EHR capture: outpatient antibiotics and UTIs are under-counted (AL2; P§3.3) · L3 survey covariates collected at enrollment, after T0 (AL4) · L4 left truncation when T0 is before enrollment (AL5; D4) · L5 site heterogeneity (AL7) · L6 death under-ascertained (AL8) · L7 no validated US EHR RCC algorithm; codes carry no histology (AL9; P§8.1) · L8 ICD-9 mapping gaps (AL10; F00 table 4) · L9 hereditary RCC syndromes without specific codes cannot be excluded (P§4) · L10 nitrofurantoin is not a pure no-drug reference (NTP TR-341; P§1.1) · L11 birth year only (O6): age at T0 can be off by up to a year, which blurs the age ≥ 40 cut-off · L12 drug capture varies by site: 309 of 1,554 EHR site-years have visits but no drug records (F01.11) · L13 some persons have more than one primary death record in v9; the earliest date is used.

## 6. Run log

| Run ID | Date (UTC) | Code hash | Commit | CDR | Cells | Ended | Output | Kept? |
|---|---|---|---|---|---|---|---|---|
| first F00+F01 | 2026-11-28 | `151a1d6eb424` | `872100a` | C2025Q4R6 | F00, F01 | ok | `workbench_outputs/2026-11-28_F00-F01_output.txt` | yes |
| F00+F01 rerun | 2026-11-28 | `a58e5386e92a` | `d368550` | C2025Q4R6 | F00, F01 | ok, but `_config`/`rutirc` were still the first run's (kernel not restarted) | `workbench_outputs/2026-11-28_F00-F01_rerun_output.txt` | partial: F00.12 and F01.11 only |
| F00+F01 final | 2026-11-28 | `2485c8c8ed67` | `c74c4dc` (PR #4) | C2025Q4R6 | F00, F01 | ok | `workbench_outputs/2026-11-28_F00-F01_final_output.txt` | yes: data checks complete |
| first blind F00–F04 | 2026-11-28 | `c4ccdaa80296` | `c793dc6` (PR #5) | C2025Q4R6 | F00–F04 | **failed** in F02 (EHR span 75.6 GiB > max_gib 60; duplicate primary death records) and F03 (nullable Int64 days supply) | none kept | no: fixed in PR #6, rerun |
| blind F00–F04 (full rerun) | 2026-11-28 | `ac8c5cb5ff51` | `PR #7` | C2025Q4R6 | F00–F04 | ok | not committed | **no:** the Q1 event total included the intermediate group while the gate did not, so the intermediate group's event count followed by subtraction (a count 1–20, by exposure group; not a primary contrast). Fixed in PR #8; F04 rerun in the same run folder. Gate (blind totals): Q1 35 events, Q2 30 events, both < 100 → D17 "descriptive + Q2b only" (O21) |

## 7. Provenance

- `pipeline/aou_core.py` is copied unchanged from BRAIN `projects/all-of-us/starter/aou_core.py` at BRAIN commit `035e4d7` (sha256 `492cdc51…`), plus the 2026-11-28 `bucket_mount` fix (no recursive walk of `~/workspace`; sha256 `7609b5c1…`), applied to BRAIN too. A fix there is a fix here too (M39).
- 2026-11-28 (PR #6), still to copy to BRAIN: `q()` normalizes BigQuery nullable Int64 / boolean columns to numpy dtypes, and `death_dates()` keeps the earliest date when v9 flags more than one primary death record per person (M40 candidate).
- Amendments after the PI signs go in `docs/AMENDMENTS-RUTI-RCC-<date>.md` (BRAIN template).
