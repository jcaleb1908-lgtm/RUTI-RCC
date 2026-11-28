> Converted to Markdown on 2026-11-27 from "Pre-cohort protocol v0.1.docx" (pandoc). Text unchanged. Status: draft for PI review, NOT signed. The signed version and its commit SHA become the prespecification stamp (STUDY-LIFECYCLE freeze line).
>
> **Revisions 2026-11-28 (team, before signature; from the first All of Us data checks, no outcome data seen):** urine culture also by SNOMED (5.3); code systems explained (Appendix A); cataract negative control limited to cataract codes, 366.11 removed (D14, A.1); nephroureterectomy kept as a confirming procedure with a sensitivity analysis without it (8.3); survey groupings for education and income and one rule for skipped answers (9.4); codes for the groups named in 9.1–9.4 written out (A.1, A.5); Appendix D sheet simplified.

**Study:** Is the kidney cancer risk seen after recurrent urinary tract
infection (UTI) explained by the infection, by the antibiotics used to
treat it, by both, or by neither (detection or reverse causation)? Does
the relationship differ in Hispanic and Puerto Rican participants?

**Data source:** NIH All of Us Research Program, Controlled Tier, data
release version 9 (C2025Q4R6; data through 1 January 2025). **Population
described:** All of Us participants with EHR data (a volunteer cohort,
not a random sample of the US; no national rates will be reported).

**Principal investigator:** Dra. Yaritza Inostroza-Nieves. **Protocol
drafting:** Valeria Román. **Workspace:** Jonathan (workspace creator).
Version 0.1, 19 November 2026. Status: draft for PI review. Written
before anyone opens the All of Us data; it contains no participant-level
data and no counts from All of Us.

## **How to read this document**

- **\[PI to decide\]** marks every open choice. Each has a proposal, the
  reason (with evidence from our literature review v3, 19 Nov 2026), and
  the alternatives. All decisions are collected in Section 16.

- **\[NOT CONFIRMED\]** marks a statement we could not check against a
  source. **\[VERIFIED\]** marks codes checked against an official file
  (sources in Appendix C).

- Analyses are labelled **primary**, **sensitivity**, **exploratory** or
  **descriptive**. Rules for choosing between definitions are fixed
  here, before any result is seen (Section 12.1).

- No All of Us concept IDs appear in this document. Each concept set is
  built inside the Researcher Workbench from the named starting concepts
  listed here, plus everything under them, and each one is checked by
  name (Appendix D). Concept IDs are never copied from a search screen.

- No analysis code. The only code is the causal diagram in dagitty
  syntax (Section 9), which contains no data.

## **Contents**

1 Research questions · 2 Target trial specification · 3 Eligibility · 4
Exclusions · 5 Recurrent UTI · 6 Antibiotics · 7 The two 2-year rules ·
8 Kidney cancer outcome · 9 Confounders and causal diagram · 10 Bias
plan · 11 Hispanic / Puerto Rican analysis · 12 Analysis plan · 13
Feasibility counts · 14 Practical checklist · 15 Team, journals,
reporting guidelines · 16 Decisions for the PI and open questions ·
Appendices A–D (code lists, references, verification log, concept-set
check sheet).

# **1. Research questions**

The literature review found no kidney-cancer study that separates the
infection from the antibiotics. Recurrent UTI has been studied once
(Mortensen 2026, Denmark: hospital-coded recurrent cystitis, kidney
cancer SIR 3.84 in year 1, falling to 1.24 (0.96–1.58) beyond 2 years,
no antibiotic or lifestyle data). Antibiotic studies show modest
associations (Park 2024 preprint, aHR 1.24 for ≥365 antibiotic days;
Boursi 2015, AOR 1.6 for \>5 penicillin courses; Petrelli 2019 pooled OR
1.28) without indication data. Neither exposure has been studied by
Hispanic ethnicity or Puerto Rican origin. Every positive study shows a
large early excess that fades, consistent with detection and protopathic
bias (Søgaard 2019; Jansåker 2025; Funston 2024).

## **1.1 Questions and the comparison that isolates each component**

| **Question**                                                                                                 | **Comparison**                                                                                                                                                                                                               | **What it isolates**                                                                                               | **What it does not isolate**                                                                                                                                                    |
|--------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Q1 — Infection.** Does recurrent UTI raise the risk of kidney cancer (RCC) compared with a single UTI?     | People with recurrent UTI vs people with a single UTI episode over the same 12-month window after a first UTI (both groups have urinary symptoms, a diagnosis and a prescriber contact).                                     | Repeated infection (number and type of episodes), with healthcare contact for urinary symptoms held roughly equal. | Recurrent UTI also brings more antibiotics, so Q1 estimates the **total** effect of recurrence (infection + its antibiotics). Q1b (below) approximates the infection-only part. |
| **Q1b — Infection with minimal systemic antibiotic** (secondary)                                             | Same as Q1, restricted to people whose UTI antibiotics in the window were only nitrofurantoin or fosfomycin.                                                                                                                 | Infection with little systemic and gut antibiotic exposure.                                                        | Nitrofurantoin itself caused kidney tubular tumours in male rats (NTP TR-341), so this is not a pure "no-drug" stratum.                                                         |
| **Q2 — Antibiotics.** For the same UTI episode, does the antibiotic class prescribed change the risk of RCC? | Active comparator, new-user design: fluoroquinolone vs nitrofurantoin/fosfomycin; TMP-SMX vs nitrofurantoin/fosfomycin; β-lactam vs nitrofurantoin/fosfomycin; plus fluoroquinolone vs β-lactam vs TMP-SMX among themselves. | The drug class, with the infection (same indication, same episode) held constant.                                  | Long-term cumulative dose; addressed by Q2b.                                                                                                                                    |
| **Q2b — Cumulative antibiotic dose** (secondary)                                                             | Cumulative courses by class over follow-up (time-varying, lagged 2 years) among people with ≥1 UTI episode, adjusted for cumulative UTI episodes.                                                                            | Dose–response for each class with infection burden adjusted.                                                       | Time-varying confounding by later infections; interpret with caution.                                                                                                           |
| **Q2c — Antibiotics without urinary infection** (exploratory)                                                | Non-urinary antibiotic courses in people with no UTI episode ever recorded.                                                                                                                                                  | Antibiotic exposure in the absence of urinary infection.                                                           | Different indications (respiratory, skin) bring their own confounding.                                                                                                          |
| **Q3 — Hispanic / Puerto Rican** (exploratory)                                                               | Q1 and Q2 estimates by group, with interaction tests (Section 11).                                                                                                                                                           | Effect modification.                                                                                               | Likely underpowered (review Section 8: a Hispanic subgroup may have ≈50–60 events **\[NOT CONFIRMED\]** estimate).                                                              |

**Interpretation key (fixed in advance).** An association in Q1 with
null Q2 contrasts points to the infection (or to detection after
infection). Q2 contrasts that are non-null within the same indication
point to the drug. Associations that disappear with the 2-year lag and
with imaging adjustment, or that also appear for the negative-control
outcome, point to detection or protopathic bias (Section 10).

> **Decision D1 — One paper or two** **\[PI to decide\]**
>
> **Proposed:** One main paper reporting Q1, Q1b, Q2 and Q2b from a
> single cohort with shared definitions; the Hispanic / Puerto Rican
> analysis (Q3) is reported inside it as exploratory. A second, short
> paper on Q3 only if the subgroup passes the event threshold in Section
> 11.4.
>
> **Reason and evidence:** The contribution identified by the review is
> the separation of infection from antibiotics (review Sections 2.3, 9);
> splitting Q1 and Q2 into two papers would separate the two halves of
> the same argument and halve the reader's ability to compare them.
> Power is limited (review Section 8.4: 300 events detect HR ≈1.7 at 10%
> exposure), so a stand-alone subgroup paper is only justified if events
> allow.
>
> **Alternatives:** (a) Two papers: Q1 (infection) and Q2 (antibiotics).
> (b) Three papers including a stand-alone Hispanic / Puerto Rican paper
> regardless of counts. (c) A methods paper on validating the RCC
> definition in All of Us (no US EHR RCC algorithm exists; review
> Section 4).

# **2. Target trial specification**

Each question is written as the randomised trial we would run if we
could, and then as its emulation in All of Us. T0 = time zero. "Outcome
counting starts" = T0 + 2 years (latency rule, Section 7) for everyone.

## **2.1 Q1 — recurrent vs single UTI (landmark design)**

| **Component**                          | **Target trial**                                                                   | **Emulation in All of Us**                                                                                                                                                                                                                                                                                                                                                                                               |
|----------------------------------------|------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Eligibility                            | Adults ≥40 with a first UTI episode, no kidney cancer, no exclusion condition.     | **Index episode (E1):** first qualifying UTI episode (Section 5) with no UTI episode in the 365 days before it. At **T0 = E1 + 365 days**: age ≥40 (completed years); ≥3 years of EHR history before T0 (Section 3.2); diagnosis and medication data available (Section 3.3); no exclusion condition in all history before T0 (Section 4); no kidney cancer code ever before T0; alive with EHR activity on or after T0. |
| Treatment strategies (exposure groups) | Not randomisable; the "trial" compares people who do vs do not develop recurrence. | **Recurrent UTI:** during \[E1, E1+365) ≥2 episodes within 180 days or ≥3 episodes within 365 days. **Single UTI (comparator):** exactly one episode (E1) during \[E1, E1+365). **Intermediate** (2 episodes \>180 days apart): excluded from the primary contrast, reported descriptively. Type: lower-only vs any upper (pyelonephritis) episode (Section 5.5).                                                        |
| Assignment                             | —                                                                                  | Status fixed at T0 using only information from the window (no future information).                                                                                                                                                                                                                                                                                                                                       |
| Time zero                              | Date of classification.                                                            | T0 = E1 + 365 days for everyone, both groups.                                                                                                                                                                                                                                                                                                                                                                            |
| Outcome                                | Incident RCC.                                                                      | Confirmed kidney cancer, ICD-10-CM C64 / ICD-9-CM 189.0 (Section 8).                                                                                                                                                                                                                                                                                                                                                     |
| Follow-up                              | From T0 until outcome, death, loss.                                                | Outcome counting starts at T0 + 2 years (people with RCC, death or last EHR activity before then are not in the analysed cohort; shown in the flow diagram). Ends at the earliest of: RCC diagnosis date, death, last EHR activity, ESRD / dialysis / kidney transplant (censoring, Section 4), 1 January 2025 data cut-off.                                                                                             |
| Causal contrast                        | Intention-to-treat analogue.                                                       | Exposure fixed at T0 (primary). Later recurrence ignored (primary); time-varying version as secondary (Section 5.6).                                                                                                                                                                                                                                                                                                     |
| Effect measure                         | Hazard ratio; risk difference.                                                     | Cause-specific HR (death competing); cumulative incidence at 5 and 10 years after outcome counting starts, and their difference.                                                                                                                                                                                                                                                                                         |
| Analysis                               | —                                                                                  | Inverse-probability-of-treatment weighted (IPTW) Cox model, weights from the Section 9 adjustment set measured before E1 (and, for UTI history, before E1); robust variance. Aalen–Johansen cumulative incidence with death as a competing event.                                                                                                                                                                        |

**Immortal time — flagged and how it is removed.** If people were
classified as "recurrent" from the date of their first UTI, the
recurrent group would be guaranteed to survive (and stay cancer-free)
until their second or third episode, which is immortal time. The
landmark design removes it: both groups must survive, without RCC, to T0
= E1 + 365 days, and classification uses only what happened inside the
window. A person who meets the rule at day 60 still starts follow-up at
day 365, like the comparator. The 2-year latency rule then moves the
start of outcome counting to E1 + 3 years for everyone.

> **Decision D2 — How Q1 handles time zero** **\[PI to decide\]**
>
> **Proposed:** Landmark at E1 + 365 days, exposure fixed at the
> landmark.
>
> **Reason and evidence:** The 12-month window is the same window as the
> recurrence definition (≥3 in 12 months; Ackerson 2024, the guideline
> definition quoted by Mortensen 2026), so classification is complete at
> the landmark and no one contributes immortal person-time. With the
> 2-year latency rule, cancers in the first 3 years after the first UTI
> are not counted, which covers the protopathic window (UTI
> prescriptions rise 9 months before renal cancer diagnosis, Funston
> 2024; largest SIRs within 3–6 months, Jansåker 2025, Søgaard 2019).
>
> **Alternatives:** (a) Sequential trials / risk-set matching: T0 = date
> the recurrence rule is met; each recurrent person matched to people
> with a UTI who have not (yet) met it, on time since E1 and calendar
> time; comparators who later become recurrent are censored with
> inverse-probability-of-censoring weights. Uses more person-time but is
> more complex. (b) Time-varying recurrent status from E1 (Cox with
> time-dependent exposure). (c) Comparator = people with no UTI, matched
> on visit count (review Section 4, secondary).

## **2.2 Q2 — antibiotic class for the same UTI episode (active-comparator, new-user design)**

| **Component**        | **Target trial**                                                                     | **Emulation in All of Us**                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|----------------------|--------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Eligibility          | Adults ≥40 with an acute UTI needing treatment; no recent use of the compared drugs. | Qualifying UTI episode (Section 5) treated with exactly one study class started within −3 to +3 days of the episode's first UTI code. At T0: age ≥40; ≥3 years of EHR history; diagnosis and medication data available; no exclusion condition; no kidney cancer ever; **no prescription of any study antibiotic class in the 90 days before T0** (new-user washout, Section 6.6). One episode per person in the primary analysis (the first eligible one). |
| Treatment strategies | Start drug A vs start drug B for this episode.                                       | Reference: nitrofurantoin or fosfomycin. Comparators: fluoroquinolone; trimethoprim ± sulfamethoxazole; β-lactam (penicillins, cephalosporins). Episodes starting ≥2 classes within the ±3-day window are excluded.                                                                                                                                                                                                                                         |
| Assignment           | Random.                                                                              | Prescriber choice; emulated by IPTW on baseline confounders, including renal function (nitrofurantoin is avoided at low eGFR), episode type and setting, prior antibiotics and prior UTIs.                                                                                                                                                                                                                                                                  |
| Time zero            | Randomisation = treatment start.                                                     | Date of the index antibiotic start.                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Outcome              | Incident RCC.                                                                        | Confirmed kidney cancer (Section 8).                                                                                                                                                                                                                                                                                                                                                                                                                        |
| Follow-up            | —                                                                                    | Outcome counting starts at T0 + 2 years; ends as in Q1.                                                                                                                                                                                                                                                                                                                                                                                                     |
| Causal contrast      | ITT.                                                                                 | ITT analogue: effect of the initial class choice, regardless of later antibiotics (later courses described, and modelled in Q2b).                                                                                                                                                                                                                                                                                                                           |
| Effect measure       | HR; risk difference.                                                                 | Cause-specific HR; cumulative incidence at 5 and 10 years.                                                                                                                                                                                                                                                                                                                                                                                                  |
| Analysis             | —                                                                                    | IPTW Cox, robust variance; balance target standardised mean difference \<0.1 for every covariate.                                                                                                                                                                                                                                                                                                                                                           |

**Immortal time in Q2:** none, provided eligibility, classification and
T0 all fall on the antibiotic start date. The ±3-day window looks back
as well as forward; a code recorded up to 3 days after the prescription
is accepted, but T0 stays at the prescription date and no outcome can
occur in that gap because outcome counting starts 2 years later.

> **Decision D3 — Reference class for Q2** **\[PI to decide\]**
>
> **Proposed:** Nitrofurantoin or fosfomycin as the reference, plus
> head-to-head fluoroquinolone vs β-lactam vs TMP-SMX.
>
> **Reason and evidence:** They are first-line for cystitis, narrow and
> urine-concentrated, so systemic and gut exposure is lowest; no earlier
> study used them as reference (review Section 4). The caveat is the NTP
> rat finding for nitrofurantoin (kidney tubular neoplasms in male
> F344/N rats, NTP TR-341), which is why the head-to-head contrasts are
> also pre-specified.
>
> **Alternatives:** (a) Non-users as reference (Boursi 2015; Park 2024)
> — mixes infection and drug. (b) Light users (1–14 days, Park 2024).
> (c) β-lactams as reference.

## **2.3 Q2b and Q2c (summary)**

- **Q2b (secondary):** cohort = everyone meeting the Q1 eligibility at
  T0 (both groups plus intermediates). Exposure = cumulative courses by
  class, updated over follow-up, counted only up to 2 years before each
  point in time (exposure lag, Section 7). Model = Cox with time-varying
  lagged exposure, adjusted for baseline confounders and lagged
  cumulative UTI episodes. A marginal structural model is the
  alternative if later infections are judged to be confounders affected
  by earlier antibiotics.

- **Q2c (exploratory):** cohort = eligible people with no UTI episode in
  all history; exposure = cumulative non-urinary systemic antibiotic
  courses (lagged 2 years); T0 = a fixed date (1 January of the first
  calendar year in which the person has ≥3 years of history and age
  ≥40).

# **3. Eligibility**

## **3.1 Age**

Age ≥40 in **completed years** at T0, computed from date of birth. If
only month and year of birth are released, day 15 is assumed and this is
stated **\[NOT CONFIRMED\]** whether the Controlled Tier gives full date
of birth.

## **3.2 "≥3 years of EHR" — measured from what**

**Definition (proposed).** EHR activity = any record from an EHR source
in the visit, condition, procedure, drug, measurement or observation
tables. Survey answers, the All of Us enrollment physical measurements,
biospecimen records and wearable data are **not** EHR activity. A person
has ≥3 years of EHR history at T0 if (a) their first EHR activity is
≥1,095 days before T0 and (b) they have ≥1 EHR visit in each of the
windows \[T0−1,095, T0−730), \[T0−730, T0−365) and \[T0−365, T0).

Reason: record span alone counts people with one visit in 2010 and the
next in 2022 as "observed"; requiring a visit in each year makes
"absence of a code" closer to "absence of the condition" for the
exclusions and confounders. The All of Us observation_period table is
not used for this rule unless we confirm how it is built **\[NOT
CONFIRMED\]**.

## **3.3 "Diagnosis and medication data available" in practice**

- **Person level:** ≥1 EHR condition record and ≥1 EHR drug record in
  the 3 years before T0.

- **Site–year level:** the person's EHR data partner contributes drug
  records in the calendar year of T0. A site–year table (drug records
  per 1,000 visits, by data partner and year) is built first; site–years
  below the 10th percentile of that ratio are flagged and excluded in a
  sensitivity analysis. Whether a data-partner identifier is available
  in v9 Controlled Tier tables is **\[NOT CONFIRMED\]**; if not, the
  person-level rule alone applies.

- **Why:** the EHR only covers care at participating health systems;
  prescriptions written elsewhere are missing, so outpatient antibiotics
  and UTIs are under-counted. Requiring active drug data reduces (does
  not remove) this misclassification.

## **3.4 Calendar window and coding eras**

Proposed T0 window: 1 January 2008 to 31 December 2022 (the last T0 that
leaves 2 years before the 1 January 2025 cut-off). ICD-9-CM was used
before 1 October 2015 and ICD-10-CM from that date, so every definition
is applied in both eras and checked in both (Section 13, Appendix D).
Coding era of T0 is a covariate and a stratification variable.

## **3.5 Time zero before vs after enrollment**

Participants enrolled from 2017 onward and had to be alive to enroll.
EHR records before enrollment are retrospective, so a T0 before
enrollment belongs only to people who survived to enroll (left
truncation). People who developed RCC and died before 2017 can never
appear.

> **Decision D4 — Allow T0 before enrollment** **\[PI to decide\]**
>
> **Proposed:** Primary: allow T0 before enrollment (retrospective EHR),
> with the survivorship issue stated. Sensitivity: delayed entry — each
> person enters the risk set at max(outcome-counting start, enrollment
> date), so only events after enrollment are counted.
>
> **Reason and evidence:** Data run to 1 January 2025 and enrollment
> began in 2017, so a prospective-only design leaves at most about 7
> years after enrollment, minus 2 years of latency; most events will be
> retrospective. Kidney cancer has comparatively good survival, so
> survivorship selection should be modest, but it is differential if
> exposure affects post-diagnosis survival. The delayed-entry analysis
> tests this (Section 12).
>
> **Alternatives:** (a) Prospective only (T0 ≥ enrollment date). (b) T0
> before enrollment but outcomes counted only after enrollment in the
> primary analysis.
>
> **Decision D5 — History rule and calendar window** **\[PI to
> decide\]**
>
> **Proposed:** ≥3 years of history with ≥1 EHR visit in each of the 3
> years before T0; T0 window 2008–2022; index UTI = first episode with
> none in the prior 365 days.
>
> **Reason and evidence:** Matches the planned "≥3 years of EHR" while
> making the look-back interpretable; about a quarter of participants
> have ≥10 years of EHR, so history is long enough for many. A 365-day
> clean window before E1 matches the recurrence window; earlier UTI
> history is adjusted for (episodes in years −3 to −1 before E1).
>
> **Alternatives:** (a) Record span only (first-to-last activity). (b)
> 3-year clean window before E1 (cleaner "new" UTI, fewer people). (c)
> T0 window starting 2010 or restricted to the ICD-10 era (T0 ≥ 1
> October 2018, so the whole 3-year look-back is ICD-10).

# **4. Exclusions**

Pre-existing disease is excluded by diagnosis codes, procedures, drugs
and labs, not by diagnosis codes alone. All codes are **\[VERIFIED\]**
against the CDC FY2027 ICD-10-CM and FY2012 (final) ICD-9-CM files, CMS
FY2027 ICD-10-PCS and CMS 2026 Physician Fee Schedule (CPT) files unless
marked; full lists in Appendix A. "History" = all EHR records before T0.
"Descendants" = all codes under a category.

| **Exclusion**                                                                  | **ICD-10-CM**                                                                                                                      | **ICD-9-CM**                                                                        | **Procedures / drugs / labs**                                                                                                                                                                                         | **Look-back** | **After T0?**                                                                                                                                                     |
|--------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Prior kidney cancer                                                            | C64.x; C65.x; Z85.52x (Z85.520, Z85.528)                                                                                           | 189.0; 189.1; V10.52; V10.53                                                        | Kidney-cancer drugs (Appendix A.4) ever before T0                                                                                                                                                                     | All history   | Cancers between T0 and T0+2y excluded by the latency rule                                                                                                         |
| Any other prior cancer except non-melanoma skin                                | C00–C96 except C44.x (Merkel cell C4A.x counts as a prior cancer); D09.0, D09.1x (urinary carcinoma in situ); Z85.x except Z85.828 | 140–208 except 173.x; 209.0–209.3x; V10.x except V10.83                             | —                                                                                                                                                                                                                     | All history   | Primary: new non-kidney cancer after T0 does **not** end follow-up. Sensitivity: censor at new cancer.                                                            |
| ESRD, dialysis, kidney transplant                                              | N18.6; Z99.2; Z49.x; Z91.15x; Z94.0; T86.1x                                                                                        | 585.6; V45.11; V45.12; V56.x **\[NOT CONFIRMED\]** all V56 sub-codes; V42.0; 996.81 | CPT 90935, 90937, 90945, 90947; 50360, 50365; PCS 0TY0\*, 0TY1\*; ICD-9 proc 39.95, 54.98, 55.69 **\[NOT CONFIRMED\]** 55.69 text. Lab: eGFR \<15 mL/min/1.73 m² on 2 tests ≥90 days apart (CKD-EPI 2021, creatinine) | All history   | **Ends follow-up** (censor) on first occurrence (acquired cystic disease and intense imaging in dialysis; immunosuppression after transplant). Sensitivity: IPCW. |
| Von Hippel–Lindau                                                              | Q85.83 (not in FY2016 file; present by FY2023). Before then VHL fell under Q85.8 (other phakomatoses).                             | 759.6 (other hamartoses, includes VHL and others)                                   | —                                                                                                                                                                                                                     | All history   | See Decision D6                                                                                                                                                   |
| Tuberous sclerosis                                                             | Q85.1                                                                                                                              | 759.5                                                                               | —                                                                                                                                                                                                                     | All history   | See D6                                                                                                                                                            |
| Polycystic kidney disease                                                      | Q61.11, Q61.19, Q61.2, Q61.3                                                                                                       | 753.12, 753.13, 753.14                                                              | —                                                                                                                                                                                                                     | All history   | See D6                                                                                                                                                            |
| Other hereditary RCC syndromes (Birt–Hogg–Dubé, hereditary leiomyomatosis RCC) | No specific code found **\[NOT CONFIRMED\]**                                                                                       | —                                                                                   | —                                                                                                                                                                                                                     | —             | Not excludable; noted as limitation                                                                                                                               |
| Chronic catheter / urinary diversion                                           | Z96.0; Z93.5x (cystostomy); Z46.6; T83.51xA (catheter-associated infection); Z97.8 not used (too broad)                            | V44.5x; V53.6; 996.64                                                               | ≥2 catheter-status or catheter-infection codes ≥90 days apart, or any cystostomy/nephrostomy status code                                                                                                              | All history   | Primary: not censored. Sensitivity: censor at first chronic-catheter code after T0.                                                                               |
| Neurogenic bladder                                                             | N31.x                                                                                                                              | 596.54; 344.61                                                                      | —                                                                                                                                                                                                                     | All history   | Not censored                                                                                                                                                      |
| Congenital urologic anomalies                                                  | Q60.x; Q62.x; Q63.x; Q64.x; Q61.4–Q61.9 (renal dysplasia, medullary cystic, other cystic) \[proposal\]                             | 753.0; 753.15–753.19; 753.2x–753.9                                                  | —                                                                                                                                                                                                                     | All history   | Not censored                                                                                                                                                      |

> **Decision D6 — Hereditary syndromes first coded after T0** **\[PI to
> decide\]**
>
> **Proposed:** Exclude only if coded before T0 (primary); sensitivity
> excludes people coded at any time.
>
> **Reason and evidence:** These conditions are present from birth, but
> they are often diagnosed because a kidney tumour was found; excluding
> people coded after T0 would remove cases because they became cases
> (conditioning on the outcome).
>
> **Alternatives:** (a) Exclude if coded at any time. (b) Exclude VHL
> ICD-9 759.6 only if a later Q85.83 code confirms it (759.6 is
> non-specific).

## **4.1 Kidney stones, immunosuppression, diabetes and CKD: exclusion, confounder or neither**

| **Condition**              | **Role**                                       | **Reason**                                                                                                                                                                                                               | **Handling**                                                                                                                                                                                                                                     |
|----------------------------|------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Kidney stones              | **Confounder** (not excluded)                  | Stones cause UTIs, cause abdominal imaging (detection), and are associated with RCC (RR 1.76, Cheungpasitporn 2015; 2.36, Bhojani 2026, partly surveillance: Hemminki 2017).                                             | Adjust (codes + stone procedures before E1 / T0). Sensitivity: exclude people with stones, as Park 2024 did.                                                                                                                                     |
| Immunosuppression          | **Confounder**                                 | Raises UTI risk and cancer risk (reduced immune surveillance). Transplant recipients are already excluded.                                                                                                               | Adjust: HIV (B20, Z21; 042), immunodeficiency due to drugs (D84.821), long-term steroids (Z79.52; V58.65), immunosuppressants (Z79.6x; V87.46 history), or ≥90 days of systemic steroid / immunosuppressant prescriptions in the year before T0. |
| Diabetes                   | **Confounder**                                 | Raises UTI risk; linked to RCC through obesity and hypertension; more healthcare contact. SGLT2 inhibitors raise genitourinary infections.                                                                               | Adjust (E08–E13; 250.x; glucose-lowering drugs; HbA1c ≥6.5%).                                                                                                                                                                                    |
| CKD (stages 1–5, not ESRD) | **Confounder before exposure; mediator after** | Baseline CKD raises UTI risk, changes antibiotic choice and is an RCC risk factor (in Hispanic people OR 2.3, Hofmann 2015). Pyelonephritis can cause scarring and CKD, which could lie on the path to RCC (review 7.1). | Adjust for CKD measured **before E1** only (N18.1–N18.5, N18.9; 585.1–585.5, 585.9; eGFR \<60 on 2 tests ≥90 days apart). Do not adjust for CKD arising after E1. ESRD excluded / censored.                                                      |

# **5. Recurrent UTI**

## **5.1 Diagnosis codes that start or extend a UTI episode**

| **Episode type**               | **ICD-10-CM (FY2027) \[VERIFIED\]**                                                                                                                                                                    | **ICD-9-CM (final) \[VERIFIED\]**                                                                                                      | **Not counted (and why)**                                                                                                                                                                                                    |
|--------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Lower UTI — cystitis**       | N30.00, N30.01 (acute cystitis ± hematuria); N30.80, N30.81 (other cystitis); N30.90, N30.91 (cystitis, unspecified)                                                                                   | 595.0 (acute); 595.89 (other: abscess of bladder, bullous, emphysematous…); 595.9 (unspecified)                                        | N30.1x interstitial, N30.2x other chronic, N30.3x trigonitis, N30.4x irradiation; 595.1, 595.2, 595.3, 595.4, 595.81, 595.82 — not acute infection                                                                           |
| **Upper UTI — pyelonephritis** | N10 (acute pyelonephritis); N12 (tubulo-interstitial nephritis, not specified as acute or chronic, where "pyelonephritis NOS" is indexed); N13.6 (pyonephrosis); N15.1 (renal and perinephric abscess) | 590.10, 590.11 (acute pyelonephritis); 590.2 (abscess); 590.80 (pyelonephritis, unspecified); 590.9 (infection of kidney, unspecified) | N11.x chronic pyelonephritis / chronic tubulo-interstitial nephritis, 590.0x chronic pyelonephritis: covariates, not episodes. 590.3 pyeloureteritis cystica, 590.81 (in diseases elsewhere, e.g., TB): not counted          |
| **Unspecified UTI**            | N39.0 (UTI, site not specified)                                                                                                                                                                        | 599.0                                                                                                                                  | See Decision D8                                                                                                                                                                                                              |
| Pregnancy / postpartum UTI     | O23.0x (kidney → upper), O23.1x (bladder → lower), O23.4x (unspecified → unspecified); O86.21 (kidney), O86.22 (bladder), O86.20/O86.29 (unspecified)                                                  | 646.6x **\[NOT CONFIRMED\]** (not checked)                                                                                             | Rare at ≥40; included for completeness                                                                                                                                                                                       |
| Not UTI episodes               | —                                                                                                                                                                                                      | —                                                                                                                                      | Z87.440 / V13.02 (history of UTI: covariate); N41.0 / 601.0 acute prostatitis (excluded; sensitivity includes it as upper-equivalent in men); R82.71 bacteriuria (asymptomatic, not counted); N34.x urethritis (not counted) |

## **5.2 Episode rule**

An episode starts on the date of a qualifying code that is ≥31 days
after the start of the previous episode. All qualifying codes from the
start date to day 30 belong to that episode (a fixed 31-day window from
the start, not a rolling window, so long chains of follow-up codes
cannot merge separate infections). The episode date is its first code
date.

> **Decision D7 — Days between codes for separate episodes** **\[PI to
> decide\]**
>
> **Proposed:** 31 days from the start of the previous episode
> (primary); 14 days (sensitivity).
>
> **Reason and evidence:** Ackerson 2024 (the US EHR recurrent-UTI
> definition) treats events within 31 days as one episode; Mortensen
> 2026 used 14 days with UTI codes or UTI antibiotics. A treatment
> failure or a follow-up visit within a month is usually the same
> infection.
>
> **Alternatives:** (a) 14 days primary. (b) Rolling window (a new
> episode only after 31 code-free days).

## **5.3 Confirmation (urine culture or antibiotic)**

A qualifying episode needs at least one qualifying code date plus either
(i) the start of an antibiotic with urinary activity (Section 6.1:
nitrofurantoin, fosfomycin, trimethoprim ± sulfamethoxazole, a
fluoroquinolone, a penicillin, a cephalosporin, a carbapenem, aztreonam,
an aminoglycoside, pivmecillinam, gepotidacin, sulopenem) within −3 to
+3 days of a code date, or (ii) a urine culture ordered or resulted
within −7 to +7 days of a code date (CPT 87086, 87088 **\[VERIFIED\]**;
or a urine-culture lab result found by name inside All of Us, coded in
LOINC or SNOMED; see Appendix A for what these code systems are).

Urine-culture lab concepts (found by name in All of Us CDR v9, 28 November
2026): **LOINC** 630-4 (Bacteria identified in Urine by Culture), 13315-7
(same, 24-hour urine), 17970-5, 17971-3, 17972-1, 17973-9, 17974-7,
44847-2, 44849-8 (Bacteria # 2 to # 8 identified in Urine by Culture);
**SNOMED** 104193001 (Bacterial culture, urine, with colony count),
104194007 (Bacterial culture, urine, with organism identification),
117010004 (Urine culture). Both systems are used because sites record the
same test in one or the other. Not used: urinalysis panels "with reflex
culture" (the culture may never have been done) and viral or
mycobacterial urine cultures.

> **Decision D8 — UTI definition details** **\[PI to decide\]**
>
> **Proposed:** (a) Code + antibiotic (±3 days) or code + urine culture
> (±7 days), adapted from Ackerson 2024. (b) Unspecified UTI (N39.0 /
> 599.0) counts towards recurrence and is classed as lower UTI unless
> the episode also contains an upper code. (c) Prostatitis excluded.
>
> **Reason and evidence:** Codes alone have PPV 84.9% for outpatient UTI
> (Germanos 2020); single codes in any position have sensitivity \>72.1%
> and PPV \<78.3% in US claims (Fortin 2022). Linking a treatment or a
> culture raises specificity, and makes Q1 and Q2 groups comparable (all
> treated or tested). Most outpatient "UTI, site not specified" is
> cystitis. Culture results without codes (Ackerson's third route) are
> not used because organism results may not be structured in All of Us
> **\[NOT CONFIRMED\]**.
>
> **Alternatives:** (a) Code alone (sensitivity analysis). (b)
> Culture-confirmed episodes only (sensitivity, if results are usable).
> (c) Require antibiotic (no culture route). (d) Put N39.0 in a separate
> third type.

## **5.4 Recurrent-UTI rule**

Recurrent UTI = ≥2 qualifying episodes within 180 days or ≥3 within 365
days (Ackerson 2024; European guideline definition quoted by Mortensen
2026). In the Q1 landmark design the rule is applied only inside \[E1,
E1+365). The number of episodes in the window (1, 2, 3, ≥4) is also
modelled as a dose term (Mahati 2026 used counts 1, 2, 3, 4, ≥5).

## **5.5 People with both cystitis and pyelonephritis**

Episode type: an episode is **upper** if any of its codes is upper;
otherwise **lower** (lower includes unspecified). Person type at T0:
**any-upper** if ≥1 upper episode in the window; **lower-only**
otherwise. Upper dominates because pyelonephritis is the plausible route
to the kidney parenchyma and the strongest trigger for imaging (review
Sections 4, 9.3); mixing N10/N11/N12 with N30/N39.0 was a weakness of
Huang 2019.

> **Decision D9 — Combined or type-specific primary exposure** **\[PI to
> decide\]**
>
> **Proposed:** Primary: recurrent UTI of any type vs single UTI.
> Secondary (pre-specified): lower-only recurrent vs single lower;
> any-upper recurrent vs single lower; single upper vs single lower.
>
> **Reason and evidence:** Power: effects after a lag are ≈1.1–1.3 and
> need ≈700–2,600 events (review Section 8.4); the review recommends
> considering combining cystitis and pyelonephritis for the primary
> test. The type-specific contrasts answer the mechanism vs detection
> question.
>
> **Alternatives:** (a) Cystitis and pyelonephritis as two co-primary
> exposures (multiplicity correction). (b) Pyelonephritis-only primary.

## **5.6 Fixed at T0 or changing over time**

Primary: fixed at T0 (landmark). Secondary: time-varying cumulative
episode count after T0, lagged 2 years (Q2b model), so later recurrences
are counted.

# **6. Antibiotics**

## **6.1 Ingredient lists by class**

All names are RxNorm ingredient names (term type IN, except where noted)
checked in NLM RxNav on 19 November 2026 **\[VERIFIED\]**. Each class
concept set = these ingredients plus all clinical drugs under them,
systemic routes only (oral, IV, IM), checked by name inside All of Us.
RxNorm identifiers (RxCUI) are listed in Appendix A.3 for checking names
only; they are not All of Us concept IDs.

| **Class (analysis group)**                     | **RxNorm ingredients**                                                                                                                                                                                                                       |
|------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **A. Nitrofurans + fosfomycin (Q2 reference)** | nitrofurantoin; fosfomycin                                                                                                                                                                                                                   |
| **B. Trimethoprim ± sulfonamide**              | trimethoprim; sulfamethoxazole (alone or combined)                                                                                                                                                                                           |
| **C. Fluoroquinolones**                        | ciprofloxacin; levofloxacin; ofloxacin; moxifloxacin; norfloxacin; gemifloxacin; delafloxacin                                                                                                                                                |
| **D1. Penicillins**                            | amoxicillin; clavulanate (with amoxicillin); ampicillin; sulbactam (with ampicillin); penicillin V; penicillin G; dicloxacillin; piperacillin; tazobactam (with piperacillin); amdinocillin pivoxil (= pivmecillinam; RxNorm uses this name) |
| **D2. Cephalosporins**                         | cephalexin; cefadroxil; cefazolin; cefuroxime; cefaclor; cefprozil; cefpodoxime; cefdinir; cefixime; cefditoren; ceftriaxone; cefotaxime; ceftazidime; cefepime; cefoxitin; cefotetan; ceftaroline; cefiderocol                              |
| D3. Carbapenems / monobactam / penems          | meropenem; ertapenem; imipenem; aztreonam; sulopenem (oral form: sulopenem etzadroxil, term type PIN, with probenecid)                                                                                                                       |
| E. Tetracyclines                               | doxycycline; minocycline; tetracycline                                                                                                                                                                                                       |
| F. Macrolides                                  | azithromycin; clarithromycin; erythromycin                                                                                                                                                                                                   |
| G. Other systemic                              | clindamycin; metronidazole; tinidazole; gentamicin; tobramycin; amikacin; plazomicin; vancomycin (IV only); linezolid; gepotidacin                                                                                                           |
| Not counted as antibiotic courses              | methenamine (urinary antiseptic; prophylaxis flag only); non-absorbed gut agents: oral vancomycin, fidaxomicin, rifaximin, oral neomycin (recorded, not counted); topical, ophthalmic, otic, vaginal and inhaled forms of any ingredient     |
| Negative-control exposures (Section 10)        | fluconazole; terbinafine; acyclovir; valacyclovir; oseltamivir                                                                                                                                                                               |

**Where other antibiotics go:** any systemic antibacterial not listed is
counted in **total** courses and in class G ("other systemic") after
review by name; it is not analysed by class. β-lactam analyses use D1 +
D2 (+D3) combined in Q2 and separately (D1, D2) in Q2b, because Boursi
2015 found the penicillin association.

## **6.2 Cumulative measure**

**Course (primary unit).** A course = a start of a class. Records of the
same class whose start dates are ≤7 days after the end of the previous
record (end = start + days supply, or start + 7 days when days supply is
missing) are merged into one course. All inpatient administrations of a
class during one hospital stay = one course. **Cumulative days supplied
(secondary)** where days supply is recorded; completeness of days supply
in All of Us drug records is **\[NOT CONFIRMED\]** and is a feasibility
item (Section 13).

## **6.3 Dose categories**

| **Measure**                 | **Categories**                                        | **Reason**                                                                                       |
|-----------------------------|-------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| Total courses (all classes) | 0, 1, 2–5, ≥6; and continuous (per course, for trend) | Matches Boursi 2015 (0, 1, 2–5, \>5) and Kilkkinen 2008 (0–1, 2–5, ≥6).                          |
| Courses by class            | 0, 1, ≥2; and continuous                              | Class-specific counts will be sparse.                                                            |
| Cumulative days (secondary) | 0, 1–14, 15–59, ≥60                                   | Lower cut-points of Park 2024 (1–14, 15–59, 60–179…); ≥365 days is unlikely outside prophylaxis. |

## **6.4 Indication rule**

A course is **urinary-indication** if a qualifying UTI code (Section
5.1) is recorded within −3 to +3 days of its start (Ackerson 2024
window); otherwise **non-urinary**. Sensitivity: ±7 days.

## **6.5 Long-term prophylaxis, setting and the exposure window**

- **Prophylaxis** = ≥90 consecutive days covered by one urinary agent
  (nitrofurantoin, trimethoprim, TMP-SMX, cephalexin, or methenamine),
  allowing gaps ≤14 days, or ≥3 starts of the same agent within 120 days
  without UTI codes at each start, or a long-term antibiotic code (Z79.2
  / V58.62) with matching prescriptions. Counted as one course per
  90-day block in course counts; flagged as a separate binary exposure;
  sensitivity excludes people with prophylaxis before T0.

- **Inpatient vs outpatient:** both counted (inpatient intravenous
  antibiotics are systemic exposure). Sensitivity: outpatient
  prescriptions only (comparable with the literature; inpatient
  administration records may vary by site **\[NOT CONFIRMED\]**).
  Setting of each course is recorded from the visit type.

- **Same exposure window for everyone:** baseline antibiotic exposure is
  counted only in \[T0 − 3 years, T0), which every eligible person has
  (Section 3.2), not over the whole (variable-length) history. After T0,
  exposure accumulates by the same rule for everyone and is lagged 2
  years.

> **Decision D10 — Antibiotic exposure measure and rules** **\[PI to
> decide\]**
>
> **Proposed:** Courses as primary unit (7-day merge rule); categories
> 0, 1, 2–5, ≥6; indication ±3 days; inpatient + outpatient; prophylaxis
> as 90-day blocks plus a flag; 3-year baseline window for everyone.
>
> **Reason and evidence:** Courses are robust to missing days supply
> (review Section 4, design table) and match the two class-specific
> studies (Boursi, Kilkkinen). ±3 days follows the only US EHR rule
> (Ackerson). The equal window avoids people with longer histories
> looking more exposed.
>
> **Alternatives:** (a) Cumulative days primary (Park 2024). (b) Defined
> daily doses (no kidney study used them). (c) Outpatient only primary.
> (d) Whole-history baseline window.

## **6.6 New-user washout for Q2**

No prescription of any study class (A–G) in the 90 days before T0. Prior
antibiotic courses in the 3-year window are covariates.

# **7. The two 2-year rules**

| **Rule**                        | **What it does**                                                                                                        | **Applies to**                                                                                                                                 | **Fixed-T0 analyses (Q1, Q2)**                                                                                                    | **Time-varying analyses (Q2b, Q2c)**                                                   |
|---------------------------------|-------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------|
| **Exposure lag (2 years)**      | Exposure that happens in the 2 years before any point in follow-up does not count at that point.                        | **Everyone's person-time**, not only cancer cases.                                                                                             | Exposure is fixed by T0, and outcomes are counted only from T0 + 2 years, so every exposure is ≥2 years old when counting starts. | At each time t, cumulative courses and episodes are those recorded before t − 2 years. |
| **Latency exclusion (2 years)** | Kidney cancers diagnosed in the first 2 years after T0 are not counted, and that person-time is not counted for anyone. | Everyone: people with RCC, death or last EHR activity before T0 + 2 years are excluded from the analysed cohort (counted in the flow diagram). | Outcome counting starts at T0 + 2 years.                                                                                          | Same, from T0.                                                                         |

**Why both are needed and why they apply to everyone:** if only cases
had their recent exposure removed, exposure would be under-counted in
cases and not in non-cases, which biases towards the null or in either
direction. Removing the same 2 years from everyone's exposure history
and everyone's follow-up keeps the comparison fair.

**Sensitivity versions:** both rules at 1 year and at 3 years (varied
together), and each varied alone. The literature supports 2 years as
above every observed protopathic window (Funston 2024: 9 months;
Mortensen 2026: kidney SIR near 1 by year 2), while the persistent late
excess in women (Jansåker 2025: 1.38 beyond 5 years) and in hematuria
(Nørgaard 2018) is why a longer window is worth testing.

> **Decision D11 — Lag and latency lengths** **\[PI to decide\]**
>
> **Proposed:** 2 years primary for both rules; 1-year and 3-year
> sensitivity versions; women reported separately.
>
> **Reason and evidence:** As above; Boursi 2015 used 1 year, Pottegård
> 2020 used 2 years, Park 2024 used washouts of 1, 3 and 5 years.
>
> **Alternatives:** Add a 5-year version (the review proposed 1 and 5
> years; it costs follow-up).

# **8. Kidney cancer outcome**

## **8.1 Codes**

- **RCC (primary outcome):** ICD-10-CM C64.1, C64.2, C64.9 (malignant
  neoplasm of kidney, except renal pelvis); ICD-9-CM 189.0 (kidney,
  except pelvis) **\[VERIFIED\]**. Codes do not carry histology; about
  nine in ten kidney (C64) cancers are RCC in registry data **\[NOT
  CONFIRMED\]** for All of Us.

- **Not RCC:** C7A.093 malignant carcinoid of kidney (ICD-9 209.24) —
  not an outcome; D30.0x (benign) and D41.0x (uncertain behaviour) and
  R93.42x (abnormal imaging of kidney) — not outcomes, but used to date
  detection (Section 8.4).

- **Secondary outcomes (urothelial):** upper-tract urothelial cancer
  C65.x, C66.x (189.1, 189.2); bladder cancer C67.x (188.x). Analysed
  separately; they do not end RCC follow-up.

## **8.2 Confirmation rule and levels**

| **Level**     | **Rule**                                                                                                                                                                                                                                                                      | **Use**                                                                                                  |
|---------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| **Confirmed** | ≥2 RCC codes on different days, the second 30–365 days after the first; **or** 1 RCC code plus a confirming kidney-cancer procedure (8.3) within −90 to +180 days of that code; **or** 1 RCC code plus a kidney-cancer systemic drug (Appendix A.4) within 365 days after it. | Primary outcome                                                                                          |
| **Possible**  | Exactly one RCC code, no confirmation.                                                                                                                                                                                                                                        | Sensitivity: confirmed + possible                                                                        |
| **Ambiguous** | RCC code and renal pelvis code (C65 / 189.1) within 90 days of each other with no confirming kidney-parenchyma procedure.                                                                                                                                                     | Excluded from confirmed; included in "possible"; counted in the secondary urothelial outcome sensitivity |

Reason for 30–365 days: codes on the same day or week are often
"rule-out" codes during work-up; a second code a month or more later
reflects ongoing care. The review's draft rule was ≥2 codes ≥30 days
apart or 1 code + nephrectomy/ablation; there is no validated US EHR
algorithm for incident RCC (review Section 4), so the rule is tested
against procedures and survey self-report (Section 12, validation).

## **8.3 Confirming procedures and keeping non-cancer nephrectomy out**

| **Procedure**                                   | **CPT (CMS 2026 PFS file) \[VERIFIED\]**                                                    | **ICD-10-PCS (FY2027) \[VERIFIED\]**                                            | **ICD-9-CM procedure**                                |
|-------------------------------------------------|---------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|-------------------------------------------------------|
| Radical / simple nephrectomy                    | 50220, 50225, 50230, 50234, 50236; laparoscopic 50545, 50546, 50548                         | 0TT00ZZ, 0TT04ZZ, 0TT04ZG, 0TT10ZZ, 0TT14ZZ, 0TT14ZG, 0TT20ZZ, 0TT24ZZ, 0TT24ZG | 55.51, 55.52; 55.54 **\[NOT CONFIRMED\]**             |
| Partial nephrectomy                             | 50240; laparoscopic 50543                                                                   | 0TB00ZZ, 0TB03ZZ, 0TB04ZZ, 0TB10ZZ, 0TB13ZZ, 0TB14ZZ (excision, non-diagnostic) | 55.4                                                  |
| Ablation of renal mass                          | 50250 (open cryo), 50542 (laparoscopic), 50592 (percutaneous RF), 50593 (percutaneous cryo) | 0T500ZZ, 0T503ZZ, 0T504ZZ, 0T510ZZ, 0T513ZZ, 0T514ZZ                            | 55.32–55.35, 55.39 **\[NOT CONFIRMED\]** exact titles |
| Renal biopsy (supports, does not confirm alone) | 50200 (percutaneous), 50205 (open)                                                          | 0TB03ZX, 0TB13ZX (diagnostic excision)                                          | 55.23, 55.24                                          |

**Nephroureterectomy** (removal of the kidney and the whole ureter: CPT
50234, 50236, 50548; ICD-9-CM 55.51) is done mostly for renal pelvis or
ureter cancer. It stays a confirming procedure in the primary analysis,
because it can only confirm a person who already has an RCC (kidney)
code. **Sensitivity:** the same cohort without nephroureterectomy as a
confirming procedure; the blind feasibility run reports how many confirmed
events depend on it (masked when 1–20). ICD-10-PCS has no single
nephroureterectomy code, so this sensitivity covers CPT and ICD-9 only.

**Rule:** procedures never define a case on their own; they only confirm
an RCC code. A procedure is not counted as confirming when, within ±90
days, the person has: a kidney-donor code (Z52.4; V59.4) or donor
nephrectomy (CPT 50300, 50320, 50547); a transplant or recipient
nephrectomy (50340, 50360, 50365; PCS codes with qualifiers Z0/Z1/Z2 =
allogeneic/syngeneic/zooplastic); kidney trauma (S37.0x; 866.x); or a
stone or chronic-infection diagnosis (N20.x, N11.x, N13.6, N15.1; 592.x,
590.0x, 590.2) **without** an RCC code in the window. Pathology codes
(88305, 88307, 88309) cannot distinguish malignant from benign and are
not used.

## **8.4 Diagnosis date**

RCC diagnosis date = earliest of: first RCC code; confirming procedure;
renal biopsy that is followed by an RCC code within 90 days. **Detection
date** (sensitivity, for detection-bias analyses) = first D41.0x, D30.0x
or R93.42x code in the 180 days before the diagnosis date, if any.

## **8.5 Renal pelvis and bladder cancer**

Renal pelvis (C65) is urothelial, not RCC, and is excluded from the
primary outcome (as in Park 2024, unlike Jansåker 2025, which used ICD-7
180 including renal pelvis). Urothelial cancers of the renal pelvis,
ureter and bladder are secondary outcomes, because chronic urinary
infection has a clearer link to urothelial than to renal-parenchymal
cancer (Lien 2013; Chagneau 2021). Prior bladder or urothelial cancer is
already an exclusion (prior cancer). A new urothelial cancer after T0
does not end RCC follow-up.

> **Decision D12 — Outcome definition** **\[PI to decide\]**
>
> **Proposed:** Confirmed level (2 codes 30–365 days apart, or code +
> procedure, or code + RCC drug) as primary; C65 excluded; urothelial
> cancers as secondary outcomes.
>
> **Reason and evidence:** Review Section 4 (Park: C64 + certificate;
> Huang: ≥2 outpatient or 1 inpatient code; Boursi: ≥1 code). In All of
> Us, EHR and survey kidney-cancer cases overlap only partly (313 in
> both of 1,266 EHR and 487 survey cases; Aschebrook-Kilfoy 2022), so a
> stricter EHR rule is preferred and validated against procedures.
>
> **Alternatives:** (a) ≥1 code (most sensitive). (b) ≥2 codes on
> different days, any interval. (c) Include C65 (comparability with
> Nordic registries).

# **9. Confounders and the causal diagram**

The diagram is drawn before any analysis, and the main adjustment set is
read from it. Paste the code below into dagitty.net (or the R package
dagitty) to reproduce the minimal sufficient adjustment sets. Node names
ending in 0 are measured before E1 (Q1) or before T0 (Q2); names ending
in 1 occur after exposure.

## **9.1 Causal diagram (dagitty code)**

dag {

rUTI \[exposure\]

RCC_dx \[outcome\]

RCC_true \[latent\]

RCC_preclinical \[latent\]

Age -\> { rUTI RCC_true HTN0 DM0 CKD0 Visits0 }

Sex -\> { rUTI RCC_true Visits0 Imaging0 }

SES -\> { Smoking Obesity Visits0 Imaging0 rUTI }

Ethnicity -\> { SES Visits0 Imaging0 Abx_win }

Smoking -\> { RCC_true Visits0 }

Obesity -\> { RCC_true DM0 HTN0 rUTI }

HTN0 -\> RCC_true

DM0 -\> { rUTI RCC_true Visits0 }

CKD0 -\> { rUTI RCC_true Abx_win Imaging0 }

Stones0 -\> { rUTI Imaging0 RCC_true }

Immuno0 -\> { rUTI RCC_true }

PriorUTI -\> { rUTI Abx_prior }

Abx_prior -\> RCC_true

Visits0 -\> { rUTI Imaging0 RCC_dx }

Imaging0 -\> RCC_dx

Site -\> { rUTI Abx_win Imaging1 RCC_dx }

Calendar -\> { rUTI Abx_win Imaging1 RCC_dx }

rUTI -\> { Abx_win CKD1 Imaging1 RCC_true }

Abx_win -\> RCC_true

CKD1 -\> RCC_true

Imaging1 -\> RCC_dx

RCC_true -\> RCC_dx

RCC_preclinical -\> { rUTI RCC_dx } // protopathic path: blocked by
design (2-year lag + latency)

}

**For Q2** set exposure = Abx_class (the class chosen for the index
episode), add nodes EpisodeType (upper vs lower), EpisodeSetting
(outpatient / emergency / inpatient), Allergy0 (penicillin allergy
Z88.0, which affects choice only) and ResistantOrg0 (prior resistant
organism), with EpisodeType, EpisodeSetting, CKD0, PriorUTI, Abx_prior,
DM0, Immuno0, Age and Sex → Abx_class.

## **9.2 Main adjustment set (from the diagram)**

**Q1:** age, sex, calendar year of T0, EHR site (stratum), smoking, BMI,
hypertension, diabetes, baseline CKD (before E1), kidney stones,
immunosuppression, number of visits in the year before E1, number of
abdominal imaging studies in the 3 years before E1, prior UTI episodes
(years −3 to −1 before E1), prior antibiotic courses (3 years before
E1), education and income (survey), self-reported race and ethnicity.
**Q2:** the same, measured before T0, plus episode type and setting and
prior resistant organism.

**Not adjusted for (mediators or descendants of exposure):** antibiotics
during the Q1 window (Abx_win: part of the total effect of recurrence;
adjusted only in Q1b by restriction), CKD arising after E1 (CKD1),
imaging after T0 (Imaging1: on the detection path; used only in the
detection-bias sensitivity analysis, Section 10), BMI after T0
(adiposity may mediate antibiotic effects: 4.3–13.2%, Liu 2026).

## **9.3 Sex**

Sex assigned at birth from The Basics survey; if missing, EHR-recorded
sex. Primary models include both sexes, adjusted for sex, with
pre-specified sex-stratified estimates and a sex × exposure interaction
test, because recurrent UTI in men is "complicated" by definition, the
antibiotic association was stronger in women in Korea (women 1.58 vs men
1.15; P-interaction 0.026, Park 2024), and the late excess after
recurrent cystitis was in women (≥2 years SIR 1.39, Mortensen 2026).

> **Decision D13 — Sex and race/ethnicity in the main models** **\[PI to
> decide\]**
>
> **Proposed:** Both sexes in the primary analysis (adjusted, plus
> stratified estimates). Self-reported race/ethnicity included in the
> adjustment set as a marker of structural factors that affect
> prescribing, imaging and risk factors.
>
> **Reason and evidence:** See 9.3 for sex. Hispanic patients receive
> more antibiotics per ambulatory visit (Young 2023) and access to
> imaging differs (Deng 2026), so race/ethnicity can confound
> exposure–detection relations.
>
> **Alternatives:** (a) Women-only primary (the recurrent-UTI definition
> comes from women; Ackerson 2024). (b) Race/ethnicity not adjusted
> (only stratified).

## **9.4 Where each variable comes from**

| **Variable**                       | **Source**                                                                                                                                                  | **Definition (window relative to T0, or E1 for Q1)**                                |
|------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|
| Age, sex                           | EHR person record; The Basics survey (sex at birth)                                                                                                         | At T0                                                                               |
| Race, ethnicity, Hispanic subgroup | **Survey only** (The Basics); never EHR race fields                                                                                                         | Answered once at or after enrollment                                                |
| Education, income, insurance       | Survey (The Basics)                                                                                                                                         | Education: less than high school / high school or GED / some college / college graduate / advanced degree. Income: <$25k / $25–50k / $50–100k / $100–200k / ≥$200k. Missing if not answered; multiple imputation (Section 12)                                |
| Smoking                            | Survey (Lifestyle) + EHR (F17.2x, Z72.0, Z87.891; 305.1, V15.82)                                                                                            | Ever vs never; EHR codes in all history before T0; survey answer regardless of date |
| BMI                                | EHR measurement (closest to T0 in the 3 years before); Z68.3x/Z68.4x, V85.3x/V85.4x as fallback; enrollment physical measurement only if within the window  | Categories \<25, 25–29.9, 30–34.9, ≥35                                              |
| Hypertension                       | EHR: I10–I16; 401–405; antihypertensive drug                                                                                                                | Any time before T0 (codes) or ≥2 prescriptions in the 3 years before T0             |
| Diabetes                           | EHR: E08–E13; 250.x; glucose-lowering drug; HbA1c ≥6.5%                                                                                                     | Any time before T0                                                                  |
| Baseline CKD                       | EHR: codes and eGFR (creatinine, CKD-EPI 2021)                                                                                                              | Before E1 (Q1) / T0 (Q2)                                                            |
| Kidney stones                      | EHR: N20.x, Z87.442; 592.x, V13.01; CPT 50080, 50081, 50590, 52352, 52353, 52356                                                                            | Any time before E1 / T0                                                             |
| Immunosuppression                  | EHR codes and drugs (Section 4.1)                                                                                                                           | Year before T0 (drugs); any time (HIV)                                              |
| Visits in the prior year           | EHR visit records                                                                                                                                           | Distinct visit days in \[T0−365, T0) (Q1: before E1)                                |
| Abdominal imaging                  | EHR procedures: CPT 74150, 74160, 74170, 74176, 74177, 74178, 74181, 74182, 74183, 76700, 76705, 76770, 76775, 72192, 72193, 72194; ICD-9 proc 87.71, 88.01 | Count in the 3 years before E1 / T0 (and time-varying after T0 for sensitivity)     |
| Prior UTI, prior antibiotics       | EHR (Sections 5, 6)                                                                                                                                         | 3 years before E1 / T0                                                              |
| Site, calendar year, coding era    | EHR metadata                                                                                                                                                | At T0                                                                               |

**Survey answers, every question:** "Skip", "Prefer not to answer" and
"Don't know" are one level, **"Not answered"**. It is never merged into
another answer.

# **10. Bias plan**

| **Bias**                            | **Why it matters here**                                                                                                                                                                                                                                                                   | **Design response**                                                                                                         | **Analysis**                                                                                                                                                                   |
|-------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Detection (surveillance)**        | More than half of RCCs are found incidentally on imaging (Rossi 2018). UTI, especially pyelonephritis, and stones lead to imaging (Hemminki 2017). In All of Us itself, UTI and other diagnoses were followed by excess diagnoses of an unrelated outcome, highest in year 1 (Wang 2025). | Single-UTI comparator (Q1); same indication (Q2); 2-year latency.                                                           | Adjust for baseline visit count and abdominal imaging; sensitivity adding time-varying post-T0 imaging; imaging-detected negative-control outcome; stratify by pyelonephritis. |
| **Protopathic (reverse causation)** | Occult RCC can present as hematuria or UTI-like symptoms; UTI prescriptions rise ≈9 months before renal cancer diagnosis (Funston 2024).                                                                                                                                                  | 2-year lag and latency (1- and 3-year versions).                                                                            | Sensitivity excluding people with hematuria (R31.x, N02.x; 599.7x) in the 2 years before T0 or with N30.01/N30.91 episodes.                                                    |
| **Confounding by indication**       | Antibiotics are given because of infection; infection-related comorbidity travels with them.                                                                                                                                                                                              | Active comparators within the same indication (Q2); recurrent vs single UTI (Q1); indication rule (urinary vs non-urinary). | IPTW with balance checks; negative-control exposure.                                                                                                                           |
| **Negative-control exposure**       | Detects residual confounding by healthcare use and frailty.                                                                                                                                                                                                                               | Antifungals and antivirals (fluconazole, terbinafine, acyclovir, valacyclovir, oseltamivir), used by Boursi 2015.           | Cumulative courses analysed exactly like antibiotics (Q2b model); expected HR ≈1.                                                                                              |
| **Negative-control outcome**        | Detects residual detection and confounding bias.                                                                                                                                                                                                                                          | See D14.                                                                                                                    | Same models with the NCO as outcome; expected HR ≈1.                                                                                                                           |
| **Unmeasured confounding**          | Diet, occupational exposures, family history incompletely coded.                                                                                                                                                                                                                          | —                                                                                                                           | E-value for each primary HR and its confidence limit closer to 1 (rare outcome, so HR approximates RR).                                                                        |
| **Death (competing risk)**          | Deaths are under-recorded in All of Us EHR.                                                                                                                                                                                                                                               | Death ends follow-up; people who died are kept.                                                                             | Cause-specific Cox (etiologic estimand) and Aalen–Johansen cumulative incidence; Fine–Gray subdistribution HR as sensitivity.                                                  |
| **Loss of EHR contact**             | People leave participating health systems; events after that are invisible.                                                                                                                                                                                                               | Follow-up ends at last EHR activity.                                                                                        | Restriction to well-observed people (Section 12).                                                                                                                              |
| **Left truncation / survivorship**  | Only people alive in 2017+ could enroll.                                                                                                                                                                                                                                                  | D4.                                                                                                                         | Delayed-entry sensitivity.                                                                                                                                                     |

> **Decision D14 — Negative-control outcomes** **\[PI to decide\]**
>
> **Proposed:** Two: (1) benign neoplasm of the adrenal gland
> (D35.00–D35.02; 227.0) — found on abdominal imaging, no known link to
> UTI or antibiotics; (2) age-related cataract (H25.x; 366.10,
> 366.12–366.19, i.e. only codes whose ICD title is a cataract; 366.11
> pseudoexfoliation of lens capsule is not a cataract and is left out) — not
> found by abdominal imaging, needs specialist contact.
>
> **Reason and evidence:** The review proposed one imaging-detected and
> one non-imaging control and noted the imaging-detected choice is ours
> **\[NOT CONFIRMED\]** as prior use. Melanoma (Boursi 2015) was
> considered for the non-imaging control but photosensitising
> antibiotics (doxycycline, fluoroquinolones) make it a poorer negative
> control.
>
> **Alternatives:** (a) Melanoma C43 (precedent: Boursi 2015). (b)
> Abdominal aortic aneurysm (I71.40–I71.43; 441.4) as imaging-detected
> control, but it shares smoking as a cause. (c) Hip osteoarthritis
> (M16.x).

# **11. Hispanic / Puerto Rican analysis**

## **11.1 Definitions (all self-reported survey answers, except residence)**

| **Group**             | **Definition**                                                                                         | **Source**                                                                                                              |
|-----------------------|--------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| Hispanic              | Selected "Hispanic, Latino, or Spanish" in the race/ethnicity question, alone or with any other answer | The Basics survey **\[NOT CONFIRMED\]** exact answer wording in v9 codebook                                             |
| Puerto Rican          | Selected the Puerto Rican subgroup answer                                                              | The Basics survey Hispanic subgroup question **\[NOT CONFIRMED\]** that the subgroup is released in v9 Controlled Tier  |
| Born in Puerto Rico   | Birthplace answer = Puerto Rico, if offered                                                            | The Basics survey **\[NOT CONFIRMED\]** whether Puerto Rico is a separate birthplace option                             |
| Living in Puerto Rico | 3-digit ZIP of current address = 006, 007 or 009                                                       | Controlled Tier ZIP3 **\[NOT CONFIRMED\]** (USPS ZIP prefixes not checked against a USPS file; 008 = US Virgin Islands) |
| Reference             | Non-Hispanic White (self-reported White, not Hispanic)                                                 | The Basics survey                                                                                                       |

## **11.2 Comparisons**

- \(a\) Hispanic vs non-Hispanic White; (b) Puerto Rican (any residence)
  vs non-Hispanic White; (c) Puerto Ricans living on the island vs
  Puerto Ricans living on the mainland (descriptive unless counts
  allow). People who skipped the question form a separate "not answered"
  group; they are not dropped from the main analysis.

- Whether participating Puerto Rico sites contribute EHR data, and for
  which years, is a feasibility item **\[NOT CONFIRMED\]**. Island
  participants without EHR data are outside the study population by
  definition.

## **11.3 Interaction tests**

- **Multiplicative:** exposure × group product term in the IPTW Cox
  model (Wald test, α = 0.05, exploratory, no multiplicity correction
  claimed).

- **Additive:** relative excess risk due to interaction (RERI) from the
  Cox model with a common reference (non-Hispanic White, unexposed),
  with delta-method or bootstrap confidence intervals; plus differences
  in 10-year cumulative incidence (exposed − unexposed) by group.

## **11.4 Minimum number of cases**

> **Decision D15 — When the subgroup becomes descriptive only** **\[PI
> to decide\]**
>
> **Proposed:** Group-specific models are fitted only if the group has
> ≥100 confirmed RCC events after all filters **and** ≥21 events in each
> compared exposure category; otherwise the group is reported
> descriptively (rates per 1,000 person-years with every cell 1–20
> suppressed).
>
> **Reason and evidence:** 100 events detect only HR ≈2.5 at 10%
> exposure (review Section 8.4), so fewer events cannot support a model;
> 21 is the smallest count the disclosure rule allows us to show. The
> reference group should be non-Hispanic White, as the largest and most
> commonly reported comparator (US Hispanic and White kidney cancer
> rates are similar, RR 0.97 men; island Puerto Rico rates are roughly
> half; review Section 6.1).
>
> **Alternatives:** (a) ≥50 events. (b) Reference = all non-Hispanic
> participants. (c) Bayesian shrinkage estimates for small groups.

# **12. Analysis plan**

## **12.1 Rules for choosing between definitions (fixed now)**

- The primary definitions are those in Sections 5–8. They are not
  changed after results are seen. A definition is changed before outcome
  analyses only if a feasibility or validation check (Section 13;
  Appendix D) shows it is broken (for example, a concept set that
  contains the wrong conditions, or a code that maps to a different
  concept in the ICD-9 era); any such change is dated and justified in a
  protocol amendment.

- Sensitivity analyses are interpreted as a group: a primary finding is
  called robust if the direction is unchanged and the HR changes by
  \<20% in the lag/latency, outcome-level and confirmation sensitivity
  analyses.

- Outcome models are run only after the cohort flow, Table 1 and balance
  diagnostics are frozen.

## **12.2 Analyses**

| **\#** | **Analysis**                                                                                                                                                                                                                                                                                    | **Label**                |
|--------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------|
| 1      | Cohort flow diagram (each step with counts; any count 1–20, or any count calculable from others, masked; complementary masking applied)                                                                                                                                                         | Descriptive              |
| 2      | Table 1 by exposure group, with a missing-data row for every variable                                                                                                                                                                                                                           | Descriptive              |
| 3      | Events, person-years and rates per 1,000 person-years by group                                                                                                                                                                                                                                  | Descriptive              |
| 4      | Cumulative incidence of RCC with death as competing risk (Aalen–Johansen), from outcome-counting start                                                                                                                                                                                          | Primary (Q1, Q2)         |
| 5      | IPTW cause-specific Cox: Q1 recurrent vs single UTI; Q2 class contrasts                                                                                                                                                                                                                         | Primary                  |
| 6      | Proportional-hazards check: scaled Schoenfeld residuals and log–log plots; if violated, period-specific HRs (2–5 and \>5 years after T0)                                                                                                                                                        | Primary (diagnostic)     |
| 7      | Q1 type-specific contrasts (D9); Q1b; Q2 head-to-head; Q2b dose–response                                                                                                                                                                                                                        | Secondary                |
| 8      | Missing data: multiple imputation by chained equations (20 imputations; outcome indicator and Nelson–Aalen cumulative hazard in the imputation model) as primary; complete-case as sensitivity                                                                                                  | Primary / Sensitivity    |
| 9      | Lag and latency at 1 and 3 years (together and separately)                                                                                                                                                                                                                                      | Sensitivity              |
| 10     | Outcome = confirmed + possible; outcome dated at detection date                                                                                                                                                                                                                                 | Sensitivity              |
| 11     | Episode window 14 days; code-only UTI; culture-confirmed UTI                                                                                                                                                                                                                                    | Sensitivity              |
| 12     | Restriction to well-observed people: ≥1 EHR visit in every year of follow-up and ≥5 years of history before T0                                                                                                                                                                                  | Sensitivity              |
| 13     | Delayed entry at enrollment (left truncation)                                                                                                                                                                                                                                                   | Sensitivity              |
| 14     | Site differences: Cox stratified by site; leave-one-site-out estimates; site-year completeness exclusion                                                                                                                                                                                        | Sensitivity              |
| 15     | ICD-10 era only (T0 ≥ 1 October 2018)                                                                                                                                                                                                                                                           | Sensitivity              |
| 16     | Excluding kidney stones; excluding hematuria; outpatient antibiotics only; censoring at new cancer; Fine–Gray                                                                                                                                                                                   | Sensitivity              |
| 17     | Detection: adding time-varying imaging after T0; stratified by pyelonephritis                                                                                                                                                                                                                   | Sensitivity              |
| 18     | Negative-control exposure and outcomes; E-values                                                                                                                                                                                                                                                | Sensitivity (bias)       |
| 19     | Sex-stratified estimates and interaction                                                                                                                                                                                                                                                        | Secondary                |
| 20     | Hispanic / Puerto Rican analyses (Section 11)                                                                                                                                                                                                                                                   | Exploratory              |
| 21     | Q2c non-urinary antibiotics in people without UTI                                                                                                                                                                                                                                               | Exploratory              |
| 22     | Validation: confirmed EHR RCC vs survey self-report of kidney cancer (concordance, by ethnicity); proportion of confirmed cases with a kidney-cancer procedure (claims-based surgery classification agreed 93% with SEER in Miller 2009); UTI episodes vs urine culture results where available | Descriptive (validation) |

> **Decision D16 — Missing data and the main estimand** **\[PI to
> decide\]**
>
> **Proposed:** Multiple imputation (primary) under missing-at-random;
> complete case as sensitivity. Main estimand: cause-specific HR with
> Aalen–Johansen cumulative incidence.
>
> **Reason and evidence:** Survey answers are often skipped and
> requiring them changes the cohort, so the main analysis must not be
> limited to survey completers. Skipping may not be at random;
> complete-case and "missing" category analyses bracket this. The
> cause-specific HR answers the causal (etiologic) question; cumulative
> incidence gives absolute risk with death competing.
>
> **Alternatives:** (a) Missing-indicator category (simple, can bias).
> (b) Complete case primary. (c) Fine–Gray subdistribution HR primary.

# **13. Feasibility counts to obtain first**

These counts are run inside the Researcher Workbench before the cohort
is built. They are for the team's planning only; any count 1–20 stays
masked even internally in shared notes. "Data Browser" = the public All
of Us Data Browser (no login; rough participant counts per code, survey
answer or drug).

| **\#** | **Count**                                                                                                                                                                                                                                                                                                                    | **Data Browser can estimate?**                                                                             |
|--------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|
| F1     | Participants with ≥1 kidney cancer code (C64 / 189.0), then at the "confirmed" level                                                                                                                                                                                                                                         | Yes for ≥1 code (participants with the condition); no for the confirmed rule                               |
| F2     | F1 after each filter in order: EHR present → age ≥40 → ≥3 years of history → exclusions → latency                                                                                                                                                                                                                            | No (needs dates and combinations)                                                                          |
| F3     | F1 and F2 among Hispanic participants and among Puerto Rican participants                                                                                                                                                                                                                                                    | Partly: survey answer counts (Hispanic; Puerto Rican subgroup if shown) but not crossed with the condition |
| F4     | People with ≥1 qualifying UTI episode; people meeting the recurrent-UTI rule; by lower-only vs any-upper                                                                                                                                                                                                                     | Partly: participants with each UTI code (N30.x, N39.0, N10); no episode or recurrence rule                 |
| F5     | Antibiotic prescriptions and courses by class; share with a UTI code within ±3 days                                                                                                                                                                                                                                          | Partly: participants with each ingredient; no courses or indication                                        |
| F6     | F1–F5 by calendar year and by coding era (ICD-9 before 1 October 2015, ICD-10 after)                                                                                                                                                                                                                                         | No                                                                                                         |
| F7     | Data-quality checks: days-supply completeness; drug record types (prescription vs inpatient administration); site–year drug completeness; data-partner identifier availability; death record source and completeness; structure of urine culture results; date-of-birth precision; EHR coverage of Puerto Rico sites by year | No                                                                                                         |
| F8     | Concept-set checks in both eras: for every definition, the standard concepts reached from the ICD-9 and ICD-10 source codes, listed by name (Appendix D)                                                                                                                                                                     | No                                                                                                         |

> **Decision D17 — Go / no-go rule after feasibility** **\[PI to
> decide\]**
>
> **Proposed:** Build the full cohort if the Q1 analysed cohort (after
> latency) is expected to have ≥300 confirmed RCC events; if 100–299,
> keep Q1 and Q2 but state the minimum detectable HR in advance; if
> \<100, report the study as descriptive plus Q2b only.
>
> **Reason and evidence:** Review Section 8.4: 300 events detect HR
> ≈1.71 (10% exposed) or 1.38 (50%); 100 events detect only HR ≈2.54 /
> 1.75. The review's illustration (≈300 events) rests on assumed numbers
> **\[NOT CONFIRMED\]**.
>
> **Alternatives:** (a) No formal rule. (b) Combine lower and upper UTI
> or drop Q1b if events are low.

# **14. Practical checklist**

| **Item**              | **Requirement / fact**                                                                                                                                                                                                          | **Status**                                              |
|-----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Data release          | Version 9, Controlled Tier C2025Q4R6, data through 1 January 2025, EHR for more than 481,000 participants. A new workspace uses the latest release; record the release in every output.                                         | Confirm the workspace shows C2025Q4R6                   |
| Tier                  | Controlled Tier needed (real dates; 3-digit ZIP; Puerto Rico is suppressed as a state in the Registered Tier).                                                                                                                  | Every analyst needs Controlled Tier access and training |
| Funding               | All workspace costs are billed to the workspace creator. Jonathan's free credits expire 24 December 2026.                                                                                                                       | See D18                                                 |
| IRB                   | Use of the Researcher Workbench is under the program's data use agreement; whether SJBSM also requires a determination is **\[NOT CONFIRMED\]**.                                                                                | See D19                                                 |
| Disclosure rule       | No count from 1 to 20 may appear, or be calculable from other numbers, in any table, figure or text (including flow diagrams, percentages, rates and person-years with small denominators, and Kaplan–Meier risk tables).       | Second person checks every output before export         |
| Data location and AI  | Participant-level data never leave All of Us; no AI tool may access them. Analysis code runs only inside the Workbench. Only aggregate results that passed the disclosure check leave it.                                       | Applies to every team member                            |
| Publication rules     | Notify All of Us ≥2 weeks before any presentation or paper; on acceptance deposit the manuscript in PubMed Central; include the program's participant acknowledgment (copy the current wording from the program at submission). | Add to project calendar                                 |
| Protocol registration | Register this protocol (e.g., OSF) before outcome analyses, so the "fixed before results" rule is visible. \[proposal\]                                                                                                         | PI                                                      |

> **Decision D18 — Who funds the workspace** **\[PI to decide\]**
>
> **Proposed:** Before 24 December 2026, move the analysis to a
> workspace billed to an institutional or grant billing account held by
> the PI (or a named funded team member), and do not build cohorts on
> the expiring credits.
>
> **Reason and evidence:** Costs follow the creator; a workspace that
> loses funding stops working mid-analysis.
>
> **Alternatives:** (a) Jonathan continues with a personal billing
> account (not recommended). (b) Apply for program credits or a research
> grant to cover compute **\[NOT CONFIRMED\]** availability.
>
> **Decision D19 — IRB** **\[PI to decide\]**
>
> **Proposed:** Ask the SJBSM IRB in writing whether a determination
> (e.g., not human subjects research) is required for analysis of
> de-identified All of Us data, and keep the answer with the protocol.
>
> **Reason and evidence:** Journals may ask for it; the institution's
> rules are not known to us.
>
> **Alternatives:** Rely on the program's agreements only.

# **15. Team, target journals and reporting guidelines**

## **15.1 Team**

| **Role**                     | **Person**                                              | **Responsibilities**                                                             |
|------------------------------|---------------------------------------------------------|----------------------------------------------------------------------------------|
| Principal investigator       | Dra. Yaritza Inostroza-Nieves                           | Decisions in Section 16; funding; IRB; final approval                            |
| Protocol and definitions     | Valeria Román                                           | Concept-set definitions and name checks; feasibility counts; manuscript drafting |
| Workspace and code           | Jonathan                                                | Workbench code, versioned; disclosure checks                                     |
| Biostatistics / epidemiology | To recruit                                              | IPTW, competing risks, multiple imputation, time-varying models                  |
| Clinical adjudication        | To recruit (urology or nephrology; infectious diseases) | Outcome and procedure rules; UTI episode and antibiotic classes                  |
| Independent checker          | To recruit                                              | Re-runs concept-set checks and disclosure checks without seeing results          |

## **15.2 Target journals (suggestions, not checked for scope or policy)**

- Cancer Epidemiology, Biomarkers & Prevention; British Journal of
  Cancer (published Mortensen 2026); Pharmacoepidemiology and Drug
  Safety (published Fortin 2022 and Kruithof 2026); JNCI Cancer
  Spectrum; Cancer Causes & Control.

- For a Puerto Rico–focused descriptive paper: Puerto Rico Health
  Sciences Journal.

## **15.3 Reporting guidelines**

- **STROBE** (cohort studies).

- **RECORD-PE** — RECORD for pharmacoepidemiology using routinely
  collected data (Langan et al., BMJ 2018;363:k3532).

- **TARGET** — Transparent Reporting of Observational Studies Emulating
  a Target Trial (Cashin et al., BMJ 2025;390:e087179; JAMA
  2025;334(12):1084). The target-trial tables in Section 2 follow its
  components.

# **16. Decisions for the PI and open questions**

## **16.1 Numbered decisions (each detailed in the section shown)**

1.  **D1 (§1)** One main paper with Q1, Q2 and exploratory Q3; a
    separate Q3 paper only if the subgroup passes D15.

2.  **D2 (§2.1)** Q1 time zero: landmark at first UTI + 365 days.

3.  **D3 (§2.2)** Q2 reference class: nitrofurantoin/fosfomycin, plus
    head-to-head contrasts.

4.  **D4 (§3.5)** Allow time zero before enrollment (primary), delayed
    entry as sensitivity.

5.  **D5 (§3)** History rule (≥1 EHR visit in each of 3 prior years), T0
    window 2008–2022, 365-day clean window before the index UTI.

6.  **D6 (§4)** Hereditary syndromes excluded only if coded before T0.

7.  **D7 (§5.2)** 31-day episode window (14 days sensitivity).

8.  **D8 (§5.3)** UTI = code + antibiotic (±3 d) or + urine culture (±7
    d); N39.0 counted as lower UTI; prostatitis excluded.

9.  **D9 (§5.5)** Combined recurrent UTI as primary exposure;
    cystitis/pyelonephritis contrasts secondary.

10. **D10 (§6)** Courses as the antibiotic unit; categories 0, 1, 2–5,
    ≥6; inpatient + outpatient; prophylaxis rule; 3-year baseline
    window.

11. **D11 (§7)** 2-year lag and latency; 1- and 3-year sensitivity (add
    5 years?).

12. **D12 (§8)** Confirmed RCC rule; C65 excluded; urothelial cancers
    secondary.

13. **D13 (§9)** Both sexes in primary analysis; race/ethnicity in the
    adjustment set.

14. **D14 (§10)** Negative-control outcomes: benign adrenal neoplasm and
    age-related cataract.

15. **D15 (§11.4)** Subgroup models only with ≥100 confirmed events and
    ≥21 per exposure category; reference non-Hispanic White.

16. **D16 (§12)** Multiple imputation primary; cause-specific HR plus
    cumulative incidence.

17. **D17 (§13)** Go / no-go thresholds (≥300 / 100–299 / \<100 events).

18. **D18 (§14)** Workspace funding before 24 December 2026.

19. **D19 (§14)** SJBSM IRB determination.

20. **D20 (§4)** ESRD, dialysis or kidney transplant after T0 ends
    follow-up (censoring) rather than being ignored.

21. **D21 (§15)** Target journal and whether to register the protocol on
    OSF before outcome analyses.

## **16.2 Open questions (to answer from inside All of Us during feasibility)**

1.  Does the v9 Controlled Tier give full date of birth, or month and
    year only?

2.  Is an EHR data-partner (site) identifier available for each record?

3.  How complete is days supply in drug records, and can prescriptions
    be told apart from inpatient administrations?

4.  Are urine culture results (organism, colony count) structured and
    usable?

5.  Where do death records come from in v9, and how complete are they
    relative to EHR?

6.  Do participating Puerto Rico health systems contribute EHR data, and
    from which years?

7.  Exact wording and release status of the Hispanic subgroup and
    birthplace answers in The Basics (codebook).

8.  How the All of Us observation-period table is built (we do not rely
    on it until this is known).

9.  Which standard (SNOMED) concepts the ICD-9 source codes in each
    definition map to, and whether any map wrongly (Appendix D).

10. Is there a kidney-cancer stage, size or histology field in any EHR
    or survey table?

# **Appendix A. Code lists**

**Code systems used in these definitions.** **ICD-9-CM / ICD-10-CM**:
diagnosis codes from billing (ICD-9 until September 2015, ICD-10 after).
**CPT** and **ICD-10-PCS / ICD-9-CM procedures**: procedure codes
(outpatient and inpatient). **RxNorm**: the US list of drug ingredients
and products (NLM). **ATC**: the WHO drug-class codes (e.g. C09 = drugs
acting on the renin–angiotensin system), used here to collect every drug
in a class. **LOINC** (Logical Observation Identifiers Names and Codes,
Regenstrief Institute): codes for laboratory tests and clinical
measurements, e.g. 2160-0 = creatinine in serum. **SNOMED CT**: a large
clinical terminology; All of Us (OMOP) uses it as the standard vocabulary
for diagnoses, and some sites also record lab tests and procedures with
it. In All of Us every ICD code "maps to" a SNOMED concept; this protocol
names ICD, CPT, LOINC and SNOMED codes and checks the mapped concepts by
name (Appendix D).

Codes not already given in full in Sections 4–9. Category codes (e.g.,
C64) mean the category and all codes under it. Sources and the exact
files used are in Appendix C.

## **A.1 Diagnoses (ICD-10-CM FY2027 / ICD-9-CM final)**

| **Concept**                                                 | **ICD-10-CM \[VERIFIED\]**                                                              | **ICD-9-CM \[VERIFIED\] unless marked**                                                |
|-------------------------------------------------------------|-----------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------|
| RCC outcome                                                 | C64.1, C64.2, C64.9                                                                     | 189.0                                                                                  |
| Renal pelvis / ureter / bladder cancer (secondary outcomes) | C65.1, C65.2, C65.9; C66.1, C66.2, C66.9; C67.0–C67.9                                   | 189.1; 189.2; 188.0–188.9                                                              |
| History of kidney cancer                                    | Z85.520, Z85.528                                                                        | V10.52, V10.53 (renal pelvis)                                                          |
| Kidney carcinoid (not outcome)                              | C7A.093                                                                                 | 209.24                                                                                 |
| Renal mass, not yet cancer (detection date only)            | D30.00–D30.02; D41.00–D41.02; R93.421, R93.422, R93.429                                 | —                                                                                      |
| Prior other cancer (exclusion)                              | C00–C96 except C44; C4A; C7A; C7B; C80; D09.0, D09.10, D09.19; Z85 except Z85.828       | 140–208 except 173; 209.0–209.3 **\[NOT CONFIRMED\]** each sub-code; V10 except V10.83 |
| Melanoma history (counts as prior cancer)                   | Z85.820                                                                                 | V10.82                                                                                 |
| Hematuria (protopathic sensitivity)                         | R31.0, R31.1, R31.21, R31.29, R31.9; N02.x; N30.01, N30.81, N30.91                      | 599.70, 599.71, 599.72                                                                 |
| Vesicoureteral reflux (covariate)                           | N13.70, N13.71, N13.721–N13.729, N13.731–N13.739; Q62.7                                 | 593.7x **\[NOT CONFIRMED\]**                                                           |
| Chronic pyelonephritis (covariate)                          | N11.0, N11.1, N11.8, N11.9                                                              | 590.00, 590.01                                                                         |
| Long-term antibiotic use                                    | Z79.2                                                                                   | V58.62                                                                                 |
| Kidney donor / trauma (exclude procedure from confirming)   | Z52.4; S37.0x                                                                           | V59.4; 866.x                                                                           |
| Obesity / BMI codes                                         | E66.x; Z68.30–Z68.39; Z68.41–Z68.45                                                     | 278.00, 278.01; V85.30–V85.39; V85.41–V85.45                                           |
| Smoking                                                     | F17.2x; Z72.0; Z87.891                                                                  | 305.1; V15.82                                                                          |
| Diabetes; hypertension                                      | E08–E13; I10, I12.x, I13.x (I11, I15, I16 **\[NOT CONFIRMED\]** not checked)            | 250.x; 401–405 (401.9 checked)                                                         |
| CKD                                                         | N18.1, N18.2, N18.30–N18.32, N18.4, N18.5, N18.9 (N18.6 = ESRD)                         | 585.1–585.5, 585.9 (585.6 = ESRD)                                                      |
| Kidney stones                                               | N20.0, N20.1, N20.2, N20.9; Z87.442                                                     | 592.0, 592.1, 592.9; V13.01                                                            |
| Immunosuppression                                           | B20; Z21; D84.821; Z79.52; Z79.6x (Z79.6x and D84.821 not in FY2016; present by FY2023) | 042; V58.65; V87.46; 279.3                                                             |
| Negative-control outcomes                                   | D35.00–D35.02; H25.x (age-related cataract)                                             | 227.0; 366.10, 366.12–366.19 (366.11 left out: not a cataract)                                             |
| Penicillin allergy (Q2 choice only)                         | Z88.0                                                                                   | —                                           
| Resistant organism (Q2 covariate)                           | Z16.x                                                                                   | V09.x                                         |
| Kidney transplant status (stops a procedure confirming RCC) | Z94.0; T86.1x                                                                           | V42.0; 996.81                                 |                                           |

## **A.2 Procedures**

Given in Sections 4, 8.3 and 9.4 (CPT from the CMS 2026 Physician Fee
Schedule relative value file, October release; ICD-10-PCS from the CMS
FY2027 codes file; ICD-9-CM procedures from the CDC FY2012 procedure
tabular). CPT is © American Medical Association; codes were checked
against the CMS file, not the AMA CPT book. Dialysis: CPT 90935, 90937,
90945, 90947. Urine culture: 87086, 87088 (87077 identification, 87186
susceptibility). Urinalysis 81001, 81003 (not used to define episodes).
Cystoscopy 52000 (detection covariate for urothelial outcomes).

## **A.3 RxNorm ingredients (name check)**

RxCUI = RxNorm concept unique identifier from NLM RxNav, 19 November
2026; for checking names only. Concept sets are built by name inside All
of Us.

| **Ingredient (RxNorm name)**         | **RxCUI** | **Ingredient**             | **RxCUI** | **Ingredient** | **RxCUI** |
|--------------------------------------|-----------|----------------------------|-----------|----------------|-----------|
| nitrofurantoin                       | 7454      | cefuroxime                 | 2194      | tetracycline   | 10395     |
| fosfomycin                           | 4550      | cefaclor                   | 2176      | azithromycin   | 18631     |
| methenamine                          | 6832      | cefprozil                  | 19552     | clarithromycin | 21212     |
| sulfamethoxazole                     | 10180     | cefpodoxime                | 20489     | erythromycin   | 4053      |
| trimethoprim                         | 10829     | cefdinir                   | 25037     | clindamycin    | 2582      |
| ciprofloxacin                        | 2551      | cefixime                   | 25033     | metronidazole  | 6922      |
| levofloxacin                         | 82122     | cefditoren                 | 83682     | tinidazole     | 10612     |
| ofloxacin                            | 7623      | ceftriaxone                | 2193      | gentamicin     | 1596450   |
| moxifloxacin                         | 139462    | cefotaxime                 | 2186      | tobramycin     | 10627     |
| norfloxacin                          | 7517      | ceftazidime                | 2191      | amikacin       | 641       |
| gemifloxacin                         | 138099    | cefepime                   | 20481     | plazomicin     | 2049549   |
| delafloxacin                         | 1927663   | cefoxitin                  | 2189      | vancomycin     | 11124     |
| amoxicillin                          | 723       | cefotetan                  | 2187      | linezolid      | 190376    |
| clavulanate                          | 48203     | ceftaroline                | 1040005   | gepotidacin    | 2709212   |
| ampicillin                           | 733       | cefiderocol                | 2265702   | fidaxomicin    | 1111103   |
| sulbactam                            | 10167     | meropenem                  | 29561     | rifaximin      | 35619     |
| penicillin V                         | 7984      | ertapenem                  | 325642    | neomycin       | 7299      |
| penicillin G                         | 7980      | imipenem                   | 5690      | fluconazole    | 4450      |
| dicloxacillin                        | 3356      | aztreonam                  | 1272      | terbinafine    | 37801     |
| piperacillin                         | 8339      | sulopenem                  | 2717833   | acyclovir      | 281       |
| tazobactam                           | 37617     | sulopenem etzadroxil (PIN) | 2717834   | valacyclovir   | 73645     |
| amdinocillin pivoxil (pivmecillinam) | 627       | doxycycline                | 3640      | oseltamivir    | 260101    |
| cephalexin                           | 2231      | minocycline                | 6980      | probenecid     | 8698      |
| cefadroxil                           | 2177      | cefazolin                  | 2180      |                |           |

## **A.4 Kidney-cancer systemic drugs (outcome confirmation; prior-cancer exclusion)**

RxNorm ingredients **\[VERIFIED\]**: sunitinib, pazopanib, axitinib,
cabozantinib, lenvatinib, sorafenib, tivozanib, everolimus,
temsirolimus, belzutifan, pembrolizumab, nivolumab, ipilimumab,
avelumab. Most are also used for other cancers, so a drug confirms RCC
only after an RCC code (Section 8.2) and never defines a case alone;
everolimus is also a transplant immunosuppressant (transplant recipients
are already excluded).

## **A.5 Drug classes and labs for covariates (ATC and LOINC)**

| **Covariate**        | **Codes**                                                                                     |
|----------------------|-----------------------------------------------------------------------------------------------|
| Antihypertensive     | ATC C02, C03, C07, C08, C09 (every drug under them)                                           |
| Glucose-lowering     | ATC A10                                                                                       |
| Systemic steroid     | ATC H02AB                                                                                     |
| Immunosuppressant    | ATC L04                                                                                       |
| Creatinine (eGFR)    | LOINC 2160-0, 38483-4                                                                         |
| HbA1c                | LOINC 4548-4, 17856-6                                                                         |
| BMI                  | LOINC 39156-5                                                                                 |
| Urine culture        | LOINC and SNOMED codes in Section 5.3                                                         |

# **Appendix B. Key references (full list in the literature review v3)**

- Ackerson BK et al. J Infect Dis 2024;230(5):e1101–e1111 (recurrent-UTI
  EHR definition).

- Aschebrook-Kilfoy B et al. PLoS One 2022;17(9):e0272522 (cancer in All
  of Us EHR vs survey).

- Boursi B et al. Eur J Cancer 2015;51(17):2655–2664 (antibiotic courses
  by class; negative controls).

- Funston G et al. Br J Gen Pract 2024;74(740):e149–e155 (pre-diagnostic
  prescriptions).

- Jansåker F et al. BMJ Public Health 2025;3(2):e002495 (cystitis and
  urogenital cancer).

- Mortensen KK et al. Br J Cancer 2026;135(10):1666–1674 (recurrent
  cystitis and urologic cancer).

- Park YJ et al. Research Square 2024, preprint (antibiotic days and
  kidney cancer, Korea).

- Wang J et al. JAMA Netw Open 2025;8(4):e256637 (detection bias in All
  of Us).

- Cashin AG et al. TARGET statement. BMJ 2025;390:e087179.

# **Appendix C. Code verification sources and log (19 November 2026)**

| **Code system**     | **Official source and file**                                                                                                                                                      | **What was checked**                                                                                                                            |
|---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| ICD-10-CM           | CDC NCHS, FY2027 code descriptions (icd10cm-order-2027.txt in icd10cm-code-descriptions-2027.zip), effective 1 October 2026; plus FY2016 and FY2023 order files to date new codes | Every ICD-10-CM code and range in this protocol: exists, title matches, billable level                                                          |
| ICD-9-CM diagnoses  | CDC NCHS, ICD-9-CM FY2012 diagnosis tabular (Dtab12), the final edition before the freeze                                                                                         | Each code's title; codes marked **\[NOT CONFIRMED\]** could not be read reliably from the file                                                  |
| ICD-9-CM procedures | CDC NCHS, ICD-9-CM FY2012 procedure tabular (Ptab12)                                                                                                                              | 55.4, 55.51, 55.52, 55.53, 55.24, 55.23, 39.95, 54.98, 87.71, 88.01 confirmed; 55.54, 55.69 and exact ablation titles 55.32–55.39 not confirmed |
| ICD-10-PCS          | CMS, 2027 ICD-10-PCS codes file (icd10pcs_codes_2027.txt)                                                                                                                         | All listed kidney resection, excision, destruction and transplant codes                                                                         |
| CPT / HCPCS         | CMS, 2026 Physician Fee Schedule relative value file, October release (RVU26D, PPRRVU2026_Oct), short descriptors                                                                 | All listed CPT codes present with matching descriptors                                                                                          |
| RxNorm              | NLM RxNav REST API (rxcui by exact name; properties), 19 November 2026                                                                                                            | Every ingredient name returns an ingredient (IN; sulopenem etzadroxil is PIN). Pivmecillinam is listed in RxNorm as "amdinocillin pivoxil"      |

Not verified against any source: Puerto Rico ZIP3 prefixes; survey
question wording; All of Us table structures and completeness (Section
16.2). All of Us facts in Sections 3 and 14 come from the team's notes
and were not re-checked here.

# **Appendix D. Concept-set check sheet (one row per definition, filled inside the Workbench)**

Lesson from our last study: a list labelled "preeclampsia" turned out to
be hyperlipidemia and coronary aneurysm. Every concept set is therefore
built from the named starting concepts in this protocol plus everything
under them, and every included concept is read by name, in both coding
eras, before use. Codes whose name does not fit are listed by the
pipeline (F00 table 2); the team's decision on each is recorded in the
study README (open items, O#).

| **Field**                   | **What to record**                                                                                            |
|-----------------------------|---------------------------------------------------------------------------------------------------------------|
| Definition name and version | e.g., "UTI episode codes v0.1"                                                                                |
| Starting concepts           | The named codes from this protocol (ICD-10-CM and ICD-9-CM source codes; RxNorm ingredients; CPT / PCS codes) |
| Descendants included?       | Yes / no, and which ones were removed after reading them                                                      |
| Standard concepts reached   | Names (not IDs) of the standard concepts the source codes map to, separately for ICD-9 and ICD-10 eras        |
| Counts by era               | Participants per era (masked 1–20), to spot era-specific mapping failures                                     |
