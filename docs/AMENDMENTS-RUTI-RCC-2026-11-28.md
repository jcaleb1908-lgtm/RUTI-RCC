# AMENDMENTS — RUTI-RCC

Prespecification stamp: the commit that recorded the PI's approval of D1–D21 (2026-11-28, PR #5). The PI's written signature on the protocol is still to be added. Entries are appended, never edited. No counts and no participant data in this file; the gate numbers are in README O21.

| # | Date | Artifact and section changed | Was | Now | Why | Estimate of the primary contrast visible when the change was made? | Stage when filed | PI approval | Commit SHA |
|---|---|---|---|---|---|---|---|---|---|
| DEV1 | 2026-11-28 | Protocol §12.1 (blind cohort description, freeze before outcome models) and STUDY-LIFECYCLE stages 7–8 | Outcome models (F05–F10) run only after the blind F00–F04 output is reviewed, committed and named in `_config.RUN["blind_frozen_output"]` | The freeze is waived. F05–F10 run in the same run folder as the F04 that built the cohorts (`_config.RUN["blinding_waived"]`, `STAGE = "final"`). The blind F00–F04 output was reviewed but not committed: its Q1 event total exposed one non-primary group's event count by subtraction (fixed in PR #8) | PI decision after the D17 gate failed, reported by Jonathan in the project thread. The PI asked for the full analysis to run now | No. Only whole-cohort event totals had been seen. One non-primary group's event count (Q1 intermediate) could be derived from the blind output; no exposure-outcome estimate had been computed | After stage 8, before stage 9 | Dra. Yaritza Inostroza-Nieves, 2026-11-28, reported by Jonathan; written confirmation pending | this commit |
| DEV2 | 2026-11-28 | D17 (go / no-go: ≥ 300 events confirmatory, 100–299 reduced, < 100 descriptive + Q2b only) | Below 100 events: descriptive analyses and Q2b only | All prespecified models run. Every estimate is labelled exploratory. An analysis where an arm or cell has ≤ 20 events stays descriptive (O19) and is marked "DESCRIPTIVE ONLY (insufficient events)" in F10.answerable | Same PI decision. The PI asked for every analysis to run, with those lacking events kept descriptive | No (as DEV1) | After stage 8, before stage 9 | as DEV1 | this commit |

## Consequences for reporting (to carry into the manuscript)

- Every hazard ratio, cumulative incidence difference and RERI is **exploratory**. None is confirmatory under the approved D17.
- The manuscript says that the blind freeze in the protocol was waived after the feasibility gate, and why (DEV1). It also says the D17 rule was overridden (DEV2).
- Analyses marked "DESCRIPTIVE ONLY" are reported as not estimable with the available events. They are never omitted (STUDY-LIFECYCLE disposition rule (a)).
- Publication: the All of Us Publication and Presentation Policy asks for program notification at least 2 weeks before a publication or presentation (BRAIN `ALL-OF-US-CONTEXT.md` §4).
