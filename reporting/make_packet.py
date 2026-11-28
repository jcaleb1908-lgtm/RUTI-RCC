"""Build the STUDY-LIFECYCLE final deliverable set from the masked ALL_TABLES.txt (outside the Workbench).

    python reporting/make_packet.py <ALL_TABLES.txt> <out_dir>

Writes (out_dir):
  01_variable_definitions.pdf       analysis variables + code lists (from _config and the protocol)
  02_demographics_socioeconomic.pdf Table 1 for Q1 and Q2 (masked, as printed by F04)
  03_cohort_flow_Q1.png/.pdf, 03_cohort_flow_Q2.png/.pdf
  04_analysis_diagram.png/.pdf      question -> analysis -> sensitivity analyses (ANALYSIS-MAP)
  05_results_by_question.pdf        one table per question (only when F05-F10 tables are present)
  06_notebook_RUTI-RCC.ipynb        the staged pipeline, outputs cleared
  07_IRB_record.md
  08_figure_pack/                   every figure with a file name, plus CAPTIONS.md
Every number comes from the masked tables; nothing is recomputed. Counts already masked stay masked.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import landscape, letter  # noqa: E402
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle  # noqa: E402
from reportlab.lib.units import inch  # noqa: E402
from reportlab.platypus import (PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table,  # noqa: E402
                                TableStyle)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "pipeline"))
import tables as T  # noqa: E402
import _config as K  # noqa: E402

SS = getSampleStyleSheet()
SMALL = ParagraphStyle("small", parent=SS["BodyText"], fontSize=7.5, leading=9)
TINY = ParagraphStyle("tiny", parent=SS["BodyText"], fontSize=6.5, leading=8)
H1, H2, BODY = SS["Heading1"], SS["Heading2"], SS["BodyText"]
STUDY = "RUTI-RCC: recurrent UTI, its antibiotics, and kidney cancer (All of Us CDR v9 C2025Q4R6)"
INK, MUTED, ACCENT, ACCENT2 = "#1f2933", "#6b7280", "#2563eb", "#0f766e"

LABELS = {
    "age_t0": "Age at T0, years", "sex": "Sex", "t0_year": "Calendar year of T0", "smoking": "Smoking (survey)",
    "bmi_cat": "BMI, kg/m2", "hypertension": "Hypertension", "diabetes": "Diabetes", "ckd_pre_e1": "Chronic kidney disease",
    "stones": "Kidney stones", "immunosuppression": "Immunosuppression", "visits_prior_year": "Visits in the prior year",
    "imaging_prior_3y": "Abdominal imaging in the prior 3 years", "prior_uti_episodes": "Prior UTI episodes",
    "prior_abx_courses": "Prior antibiotic courses", "education": "Education (survey)", "income": "Household income (survey)",
    "race_eth": "Race and ethnicity (survey)", "episode_type": "Index episode type", "episode_setting": "Index episode setting",
    "resistant_org": "Resistant organism code", "pcn_allergy": "Penicillin allergy code", "egfr_cat": "eGFR, mL/min/1.73 m2",
    "coding_era": "Coding era at T0", "years_t0_to_enrollment": "Years from T0 to enrollment",
    "ehr_years_before_t0": "Years of EHR before T0", "hispanic": "Hispanic or Latino", "puerto_rican": "Puerto Rican",
    "pr_residence": "Residence (ZIP3)", "insurance": "Insurance (survey)", "prophylaxis_before": "Prior antibiotic prophylaxis",
    "hematuria_2y": "Hematuria code in the prior 2 years",
}
GROUP_NAMES = {"recurrent": "Recurrent UTI", "single": "Single UTI", "intermediate": "Intermediate",
               "A": "A nitrofurantoin / fosfomycin", "B": "B TMP-SMX / trimethoprim", "C": "C fluoroquinolone",
               "D": "D beta-lactam"}

# Analysis variables (protocol sections in brackets). Definitions follow the approved protocol v0.1 and D1-D21.
VARIABLES = [
    ("UTI episode", "Exposure building block", "A 31-day fixed window started by a UTI code (cystitis, pyelonephritis, unspecified UTI); qualifying when a urine culture is within +/-7 days or an antibiotic within +/-3 days of the code (D7)", "condition_occurrence, observation, measurement, procedure_occurrence (ICD-9/10); drug_exposure; measurement (urine culture LOINC/SNOMED)", "P§5.1-5.3"),
    ("E1", "Index date", "The first qualifying UTI episode with a clean prior 365 days (no coded UTI episode) whose T0 meets eligibility (O1: first eligible)", "derived", "P§2.1, O1, O2"),
    ("q1_group (Q1 exposure)", "Exposure", "Episodes in [E1, E1 + 365 d): recurrent = >= 2 episodes within 180 d or >= 3 within 365 d; single = E1 only; intermediate = 2 episodes > 180 d apart (reported separately)", "derived from episodes", "P§2.1, 5.4, D9"),
    ("T0 (Q1)", "Time zero", "E1 + 365 days (landmark, D2)", "derived", "P§2.1"),
    ("q2_class (Q2 exposure)", "Exposure", "Antibiotic class started within +/-3 d of a qualifying episode: A nitrofurantoin/fosfomycin (reference), B TMP-SMX/trimethoprim, C fluoroquinolone, D beta-lactam; new users (no study antibiotic in the prior 90 d); one class only", "drug_exposure (RxNorm ingredients, systemic routes)", "P§2.2, 6.1, 6.6"),
    ("T0 (Q2)", "Time zero", "Start date of the index antibiotic", "derived", "P§2.2"),
    ("Antibiotic course", "Exposure building block", "Same-class records merged when a start is <= 7 days after the running end; one course per inpatient stay; urinary when a UTI code is within +/-3 d; prophylaxis periods counted as one course per 90-day block", "drug_exposure, visit_occurrence", "P§6.2-6.5, D10"),
    ("Cumulative courses (Q2b)", "Time-varying exposure", "Courses by class, 0 / 1 / 2-5 / >= 6, lagged 2 years", "derived", "P§2.3, D10, D11"),
    ("Confirmed RCC (event)", "Outcome", "C64 / 189.0 plus a second code 30-365 d later, a confirming procedure (nephrectomy, partial nephrectomy, ablation, nephroureterectomy) within -90..+180 d, or a kidney-cancer drug; C65 excluded (D12)", "condition_occurrence, procedure_occurrence, drug_exposure", "P§8.1-8.3"),
    ("RCC diagnosis date", "Outcome date", "Earliest RCC code (O7)", "derived", "P§8.4"),
    ("Latency", "Design", "Outcome counting starts 2 years after T0 (D11)", "derived", "P§7"),
    ("Follow-up end", "Design", "Earliest of RCC, death, last EHR activity, ESRD (censor, D20), 2025-01-01", "death, visit_occurrence, derived ESRD", "P§7"),
    ("Age at T0", "Eligibility, covariate", "Completed years; birth year only is real in v9 (O6, L11); >= 40 required", "person", "P§3.1"),
    ("EHR history", "Eligibility", ">= 3 years of EHR with a visit in each of the 3 years before T0, and >= 1 condition and >= 1 drug record in the 3 years before T0", "visit_occurrence, *_ext src_id (EHR sites only, O5)", "P§3.2-3.3, D5"),
    ("Exclusions", "Eligibility", "Before T0: kidney cancer, other cancer except non-melanoma skin, ESRD / dialysis / transplant, VHL / tuberous sclerosis / PKD, chronic catheter or diversion, neurogenic bladder, congenital urologic anomaly", "condition_occurrence and others", "P§4"),
    ("Sex", "Covariate, modifier", "Sex at birth from the Basics survey; EHR gender when missing", "survey, person", "P§9.3"),
    ("Smoking", "Covariate", "Ever = >= 100 cigarettes in lifetime (survey); never / unknown otherwise", "survey", "P§9.4"),
    ("BMI", "Covariate (imputed)", "Last value in the 3 years before E1 (Q1) or T0 (Q2); 12-90 kg/m2 kept", "measurement LOINC 39156-5", "P§9.4"),
    ("Hypertension, diabetes", "Covariates", "A diagnosis code or a drug of the class (ATC C02-C09; A10; HbA1c >= 6.5% for diabetes) before E1 / T0", "condition_occurrence, drug_exposure, measurement", "P§9.4, O9"),
    ("CKD, eGFR", "Covariates", "CKD code or two eGFR < 60 tests >= 90 d apart; eGFR by CKD-EPI 2021 from the last creatinine (Q2 adjustment)", "measurement LOINC 2160-0, 38483-4", "P§9.4"),
    ("Stones, immunosuppression", "Covariates", "Kidney stone code; immunosuppressant (ATC L04, H02AB >= 30 d) or transplant / HIV codes before E1 / T0", "condition_occurrence, drug_exposure", "P§4.1, 9.4"),
    ("Health-care use", "Covariates", "Visits in the prior year; abdominal imaging in the prior 3 years; prior UTI episodes; prior antibiotic courses", "visit_occurrence, procedure_occurrence", "P§9.2"),
    ("Education, income", "Covariates (imputed)", "Basics survey; 5 levels each; skip / prefer not / don't know = Not answered", "survey", "P§9.4, O10"),
    ("Race and ethnicity", "Covariate", "Basics survey; Hispanic of any race first; Not answered kept as a level", "survey", "P§9.4, 11.1"),
    ("Hispanic, Puerto Rican", "Effect modifiers", "Self-reported Hispanic or Latino; Puerto Rican subgroup answer", "survey", "P§11.1"),
    ("Insurance", "Descriptive", "Insurance type, first ticked in the order Medicaid, Medicare, Military or VA, Private, Other, None", "survey", "P§9.4"),
    ("EHR site", "Stratum", "src_id of the last EHR visit before T0", "visit_occurrence_ext", "P§9.2, O18"),
    ("Index episode type and setting (Q2)", "Covariates", "Upper vs lower; outpatient / emergency / inpatient from the episode visit", "derived", "P§9.2"),
]


def styled_table(data, widths=None, font=7, header_bg="#eef2ff", zebra=True):
    t = Table(data, colWidths=widths, repeatRows=1)
    st = [("FONT", (0, 0), (-1, -1), "Helvetica", font), ("FONT", (0, 0), (-1, 0), "Helvetica-Bold", font),
          ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)), ("VALIGN", (0, 0), (-1, -1), "TOP"),
          ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#94a3b8")),
          ("TOPPADDING", (0, 0), (-1, -1), 1.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5)]
    if zebra:
        for r in range(1, len(data)):
            if r % 2 == 0:
                st.append(("BACKGROUND", (0, r), (-1, r), colors.HexColor("#f8fafc")))
    t.setStyle(TableStyle(st))
    return t


def p(text, style=SMALL):
    return Paragraph(str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"), style)


def doc(path, story, land=False):
    d = SimpleDocTemplate(str(path), pagesize=landscape(letter) if land else letter, leftMargin=0.5 * inch,
                          rightMargin=0.5 * inch, topMargin=0.5 * inch, bottomMargin=0.5 * inch, title=path.stem)

    def foot(c, dd):
        c.saveState(); c.setFont("Helvetica", 7); c.setFillColor(colors.HexColor(MUTED))
        c.drawString(0.5 * inch, 0.3 * inch, f"{STUDY} | masked aggregates only (counts 1-20 shown as <=20)")
        c.drawRightString(dd.pagesize[0] - 0.5 * inch, 0.3 * inch, f"page {dd.page}")
        c.restoreState()
    d.build(story, onFirstPage=foot, onLaterPages=foot)


# ------------------------------------------------------------------------------------------------- 01
def variable_definitions(out, code_hash):
    s = [Paragraph("Variable definitions", H1), p(STUDY, BODY),
         p(f"Pipeline {K.PIPELINE_VERSION}; {K.PROTOCOL_VERSION}; code {code_hash}.", SMALL), Spacer(1, 6),
         Paragraph("A. Analysis variables", H2),
         p("Short definitions. The protocol sections named in the last column and pipeline/_config.py are authoritative.", SMALL)]
    rows = [["Variable", "Role", "Definition", "Source (OMOP tables)", "Protocol"]]
    rows += [[p(a), p(b), p(c), p(d), p(e)] for a, b, c, d, e in VARIABLES]
    s.append(styled_table(rows, [1.2 * inch, 0.9 * inch, 3.0 * inch, 1.6 * inch, 0.7 * inch]))
    s += [PageBreak(), Paragraph("B. Diagnosis code sets", H2),
          p("Resolved by code inside the CDR (F00); '*' = all descendants. Role: E exposure, O outcome, X exclusion, C covariate, D diagnostic work-up, N negative control.", SMALL)]
    rows = [["Set", "Role", "ICD-10-CM", "ICD-9-CM"]]
    for k, v in K.DX.items():
        rows.append([p(k), p(v.get("role", "")), p(", ".join(v.get("icd10", [])), TINY), p(", ".join(v.get("icd9", [])), TINY)])
    s.append(styled_table(rows, [1.5 * inch, 0.4 * inch, 3.8 * inch, 1.7 * inch]))
    s += [Spacer(1, 8), Paragraph("C. Procedure code sets", H2)]
    rows = [["Set", "CPT / HCPCS", "ICD-10-PCS", "ICD-9 procedure"]]
    for k, v in K.PX.items():
        rows.append([p(k), p(", ".join(v.get("cpt", [])), TINY), p(", ".join(v.get("pcs", [])), TINY), p(", ".join(v.get("icd9p", [])), TINY)])
    s.append(styled_table(rows, [1.5 * inch, 2.0 * inch, 2.6 * inch, 1.3 * inch]))
    s += [Spacer(1, 8), Paragraph("D. Laboratory codes", H2)]
    rows = [["Lab", "LOINC", "SNOMED", "Plausible range"]]
    for k, v in K.LABS.items():
        rng = v.get("range_mgdl") or v.get("range") or ""
        rows.append([p(k), p(", ".join(v.get("loinc", [])), TINY), p(", ".join(v.get("snomed", [])), TINY), p(rng)])
    s.append(styled_table(rows, [1.5 * inch, 3.0 * inch, 1.6 * inch, 1.3 * inch]))
    s += [Spacer(1, 8), Paragraph("E. Antibiotic classes (RxNorm ingredients, systemic routes)", H2)]
    rows = [["Class", "Ingredients"]] + [[p(k), p(", ".join(v), TINY)] for k, v in K.ABX_CLASSES.items()]
    rows += [[p("Kidney-cancer drugs"), p(", ".join(K.RCC_DRUGS), TINY)], [p("Negative-control drugs"), p(", ".join(K.NEG_CONTROL_RX), TINY)]]
    s.append(styled_table(rows, [1.5 * inch, 5.9 * inch]))
    s += [Spacer(1, 8), Paragraph("F. Survey items", H2)]
    rows = [["Question concept code", "Variable"]] + [[p(k), p(v)] for k, v in K.SURVEY_ITEMS.items()]
    s.append(styled_table(rows, [3.0 * inch, 4.4 * inch]))
    doc(out / "01_variable_definitions.pdf", s)


# ------------------------------------------------------------------------------------------------- 02
def table1_story(tab, qn):
    df = tab["df"].copy()
    ncols = [c for c in df.columns if c.startswith("n_") and c != "n_missing_all"]
    groups = [c[2:] for c in ncols]
    smd = [c for c in df.columns if c.startswith("smd_")]
    head = ["Characteristic", "Level"] + [GROUP_NAMES.get(g, g) for g in groups] + [c.replace("smd_", "SMD ").replace("_vs_", " vs ") for c in smd]
    rows, last = [head], None
    for _, r in df.iterrows():
        var = LABELS.get(r["var"], r["var"])
        cells = []
        for g in groups:
            v = r.get(f"n_{g}", "") or r.get(f"pctl_{g}", "")
            cells.append(v)
        sm = []
        for c in smd:
            x = T.num([r[c]]).iloc[0]
            sm.append("" if pd.isna(x) else f"{x:.2f}")
        rows.append([p(var if var != last else ""), p(r["level"])] + [p(c) for c in cells] + [p(x) for x in sm])
        last = var
    n_g = len(groups) + len(smd)
    w = [1.6 * inch, 1.3 * inch] + [(7.0 * inch) / n_g] * n_g
    return [styled_table(rows, w, font=6.5),
            p(tab["note"].split("|")[-2] if "smd" in tab["note"] else "", TINY)]


def demographics(out, tabs, code_hash):
    s = [Paragraph("Demographic and socioeconomic characteristics", H1), p(STUDY, BODY),
         p("Counts are persons; continuous variables are median and quartiles (p25 / p50 / p75). SMD = standardized difference "
           "before weighting (|SMD| > 0.1 = imbalance). Cells of 1-20 persons show <=20; a few larger cells also show <=20 so "
           "that a small cell cannot be derived from a total (All of Us dissemination rule).", SMALL), Spacer(1, 6)]
    for qn, title in (("Q1", "Q1: recurrent vs single UTI (landmark at E1 + 365 d)"),
                      ("Q2", "Q2: antibiotic class for the same UTI episode (new users)")):
        key = f"F04.table1 {qn}"
        if key in tabs:
            s += [Paragraph(title, H2)] + table1_story(tabs[key], qn)
            arms = tabs.get(f"F04.arms {qn}")
            if arms is not None:
                s.append(p("Analysed persons: " + "; ".join(f"{GROUP_NAMES.get(g, g)} {n}" for g, n in zip(arms["df"]["group"], arms["df"]["n"])), SMALL))
            s.append(PageBreak())
    s.append(p(f"Code {code_hash}.", TINY))
    doc(out / "02_demographics_socioeconomic.pdf", s, land=True)


# ------------------------------------------------------------------------------------------------- 03
def flow_figure(tab, qn, path_base):
    df = tab["df"]
    steps = list(df["step"])
    n = list(df["n"])
    ex = list(df["excluded"]) if "excluded" in df else [""] * len(steps)
    h = 0.62 * len(steps) + 0.6
    fig, ax = plt.subplots(figsize=(8.5, h))
    ax.set_xlim(0, 10); ax.set_ylim(0, len(steps)); ax.axis("off")
    for i, (st, nn, e) in enumerate(zip(steps, n, ex)):
        y = len(steps) - i - 0.5
        last = i == len(steps) - 1
        ax.add_patch(plt.Rectangle((0.2, y - 0.3), 5.6, 0.6, fc="#eef2ff" if not last else "#dbeafe",
                                   ec=ACCENT, lw=0.8))
        ax.text(0.35, y, st[0].upper() + st[1:], va="center", fontsize=7.2, color=INK, wrap=True)
        ax.text(5.65, y, f"n = {nn}" if nn else "n hidden", va="center", ha="right", fontsize=7.5, color=INK, weight="bold")
        if i < len(steps) - 1:
            ax.annotate("", xy=(3.0, y - 0.32), xytext=(3.0, y - 0.68), arrowprops=dict(arrowstyle="<-", color=MUTED, lw=0.7))
        if e:
            ax.text(6.1, y, f"excluded {e}", va="center", fontsize=7, color=MUTED)
    ax.set_title(f"{qn} cohort construction (masked)", fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{path_base}.{ext}", dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------- 04
QNAMES = {"ALL": "Both cohorts:\ndescription and\nsensitivity analyses", "Q1": "Recurrent vs\nsingle UTI", "Q2": "Antibiotic class\nfor the same episode",
          "Q1b": "Infection with minimal\nantibiotic", "Q2b": "Cumulative courses\nby class", "Q3": "Hispanic and\nPuerto Rican groups",
          "Q2c": "Non-urinary antibiotics,\nno UTI (not run)"}


def analysis_diagram(path_base):
    am = pd.read_csv(sorted((ROOT / "docs").glob("ANALYSIS-MAP-*.csv"))[-1], comment="#")
    qs = am["question_id"].drop_duplicates().tolist()
    fig, ax = plt.subplots(figsize=(12.5, 0.36 * len(am) + 2.4))
    ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(-6, len(am) + 1)
    y = len(am)
    for q in qs:
        sub = am.loc[am["question_id"].eq(q)]
        qlabel = QNAMES.get(q, q)
        y0 = y
        for _, r in sub.iterrows():
            parent = isinstance(r.get("parent_id"), str) and r["parent_id"].strip() != ""
            x = 3.9 if parent else 3.2
            col = ACCENT2 if parent else ACCENT
            ax.text(x, y, f"{r['analysis_id']}  {str(r['question'])[:100]}  ({r['label']})", fontsize=6.6, va="center", color=col)
            ax.text(11.9, y, str(r.get("cell", "")).replace("_not_run", " (not run)"), fontsize=6.3, va="center", ha="right", color=MUTED)
            ax.plot([3.0, x - 0.08], [y, y], color="#cbd5e1", lw=0.6)
            y -= 1
        ax.plot([3.0, 3.0], [y + 1, y0], color="#94a3b8", lw=0.8)
        ax.text(0.1, (y0 + y + 1) / 2, f"{q}\n{qlabel}", fontsize=7.3, va="center", color=INK, weight="bold", linespacing=1.3)
        y -= 0.9
    ax.set_title("Analysis map: each question, its analyses (blue) and sensitivity analyses (green), and the cell that runs it",
                 fontsize=9, loc="left", color=INK)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{path_base}.{ext}", dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------- 05 / 08
def descriptive_rows(tabs):
    from scipy.stats import chi2
    rows = []
    for q, lab in (("Q1", "Q1 cohort (recurrent + single UTI)"), ("Q2", "Q2 cohort (new users of a study antibiotic class)")):
        t = tabs.get(f"F04.events_total {q}")
        if t is None:
            continue
        r = t["df"].iloc[0]
        n, ev, py = T.num([r["n"]]).iloc[0], T.num([r["events"]]).iloc[0], T.num([r["person_years"]]).iloc[0]
        if pd.isna(ev) or pd.isna(py):
            continue
        lo = chi2.ppf(0.025, 2 * ev) / 2 if ev > 0 else 0.0
        hi = chi2.ppf(0.975, 2 * (ev + 1)) / 2
        rows.append({"Cohort": lab, "Persons": f"{int(n):,}", "RCC events": f"{int(ev)}", "Person-years": f"{int(py):,}",
                     "Rate per 1,000 PY (95% CI)": f"{1000 * ev / py:.2f} ({1000 * lo / py:.2f} to {1000 * hi / py:.2f})",
                     "Deaths (competing)": str(r.get("deaths_competing", ""))})
    return rows


def results(out, tabs, figdir, captions, code_hash):
    allr = tabs.get("F10.all")
    if allr is None:
        return False
    R = allr["df"].copy()
    ans = tabs.get("F10.answerable")
    s = [Paragraph("Results by question", H1), p(STUDY, BODY),
         p("Every estimate is exploratory: the feasibility gate (D17) needs >= 100 events for a reduced analysis and >= 300 for a "
           "confirmatory one, and the blind freeze was waived by the PI (DEV1, DEV2). A hazard ratio is shown only when each "
           "compared arm has >= 21 events; otherwise the analysis stays descriptive.", SMALL), Spacer(1, 6)]
    desc = descriptive_rows(tabs)
    if desc:
        s += [Paragraph("Descriptive results that can be reported", H2),
              p("Whole-cohort RCC incidence after the 2-year latency (all exposure groups together), with exact Poisson 95% CI; "
                "and every cumulative incidence that met the >= 21 events rule.", SMALL)]
        s += [styled_table([list(desc[0].keys())] + [[p(v) for v in r.values()] for r in desc], font=7), Spacer(1, 6)]
        cif = tabs.get("F05.cif_Q1")
        if cif is not None and "CIF_pct" in cif["df"]:
            shown = cif["df"].loc[cif["df"]["status"].eq("shown")]
            for _, r in shown.iterrows():
                s.append(p(f"Q1 {r['group']}: cumulative incidence of RCC {float(r['years_after_counting_start']):.0f} years after "
                           f"outcome counting starts = {float(r['CIF_pct']):.2f}% (95% CI {float(r['lo']):.2f} to {float(r['hi']):.2f}), "
                           "IPTW, death competing. The recurrent-UTI curve could not be shown (<= 20 events).", SMALL))
        hr = tabs.get("F08.hisp_rates")
        if hr is not None:
            d = hr["df"].loc[~hr["df"]["group"].str.contains("living in Puerto Rico")]
            s += [Spacer(1, 4), p("Person-years and RCC events by Hispanic and Puerto Rican group (A20). Rows for Puerto Ricans living "
                                  "in Puerto Rico are left out: <= 20 such persons per exposure group (Table 1).", SMALL),
                  styled_table([list(d.columns)] + d.astype(str).values.tolist(), font=6.5)]
        s.append(PageBreak())
    if ans is not None:
        s += [Paragraph("What each analysis could answer", H2)]
        a = ans["df"]
        rows = [["Analysis", "Question", "Estimates shown / planned", "Answer"]]
        rows += [[p(f"{r['analysis_id']} {r['analysis']}"), p(r["question"]), p(f"{r['estimates_shown']} / {r['estimates_total']}"),
                  p(r["answer"])] for _, r in a.iterrows()]
        s += [styled_table(rows, [4.6 * inch, 0.7 * inch, 1.3 * inch, 2.8 * inch]), PageBreak()]
    for q in [x for x in R["question"].drop_duplicates() if x]:
        sub = R.loc[R["question"].eq(q)]
        s.append(Paragraph(f"{q}", H2))
        rate = tabs.get(f"F05.rates_{q}")
        if rate is not None:
            d = rate["df"]
            s += [p(rate["title"], SMALL), styled_table([list(d.columns)] + d.astype(str).values.tolist(), font=6.5), Spacer(1, 4)]
        cols = [c for c in ["analysis_id", "analysis", "contrast", "HR", "lo", "hi", "p", "E_value", "status"] if c in sub.columns]
        rows = [cols] + [[p(v, TINY) for v in r] for r in sub[cols].astype(str).values.tolist()]
        s += [styled_table(rows, font=6.5), PageBreak()]
    for key in ("F10.robust",):
        if key in tabs and len(tabs[key]["df"]):
            d = tabs[key]["df"]
            s += [Paragraph("Robustness rule (P§12.1)", H2), styled_table([list(d.columns)] + d.astype(str).values.tolist(), font=6.5)]
    s.append(p(f"Code {code_hash}.", TINY))
    doc(out / "05_results_by_question.pdf", s, land=True)
    # forest plot of every shown HR
    if not {"HR", "lo", "hi"} <= set(R.columns):
        return True
    R["HRn"], R["lon"], R["hin"] = T.num(R["HR"]), T.num(R["lo"]), T.num(R["hi"])
    F = R.loc[R["status"].eq("shown") & R["HRn"].notna()]
    if len(F):
        F = F.iloc[:60]
        fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(F) + 1.2))
        yy = np.arange(len(F))[::-1]
        ax.errorbar(F["HRn"], yy, xerr=[F["HRn"] - F["lon"], F["hin"] - F["HRn"]], fmt="o", color=ACCENT, ms=3.5,
                    ecolor="#94a3b8", elinewidth=0.9, capsize=0)
        ax.axvline(1, color=MUTED, lw=0.8, ls="--")
        ax.set_xscale("log"); ax.set_yticks(yy)
        ax.set_yticklabels([f"{a} {str(b)[:48]} | {c}" for a, b, c in zip(F["analysis_id"], F["analysis"], F["contrast"])], fontsize=6.3)
        ax.set_xlabel("Hazard ratio (log scale), 95% CI", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("Every hazard ratio shown (exploratory)", fontsize=9.5, loc="left")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(figdir / f"fig4_forest_all_HR.{ext}", dpi=200)
        plt.close(fig)
        captions.append(("fig4_forest_all_HR", "Hazard ratios and 95% confidence intervals for every contrast with >= 21 events "
                         "in each compared arm. Exploratory (D17 not met; DEV1, DEV2). Labels carry no counts."))
    return True


# ------------------------------------------------------------------------------------------------- 06 / 07
def notebook(out):
    import nbformat as nbf
    nb = nbf.v4.new_notebook()
    nb.cells.append(nbf.v4.new_markdown_cell(f"# {STUDY}\nStaged pipeline, outputs cleared. Run the first cell, then each stage "
                                             "in order (or `RUN_IN_WORKBENCH.py` as one cell). Participant-level data never "
                                             "leave the Workbench; only masked `ALL_TABLES.txt` is copied out."))
    nb.cells.append(nbf.v4.new_code_cell((ROOT / "pipeline" / "RUN_IN_WORKBENCH.py").read_text()))
    for f in ["F00_setup", "F01_inventory", "F02_extract", "F03_phenotypes", "F04_cohorts", "F05_main", "F06_timevarying",
              "F07_missing", "F08_modification", "F09_sensitivity", "F10_results"]:
        nb.cells.append(nbf.v4.new_markdown_cell(f"## {f}"))
        nb.cells.append(nbf.v4.new_code_cell((ROOT / "pipeline" / f"{f}.py").read_text()))
    nbf.write(nb, str(out / "06_notebook_RUTI-RCC.ipynb"))


def irb(out):
    (out / "07_IRB_record.md").write_text(
        "# IRB record: RUTI-RCC\n\n"
        "| Field | Value |\n|---|---|\n"
        "| IRB of record | All of Us Research Program central IRB covers Researcher Workbench use; local determination: to be confirmed by Jonathan |\n"
        "| Local IRB number | PENDING (asked 2026-11-28) |\n| Determination / approval date | PENDING |\n"
        "| Data tier | Controlled Tier, CDR v9 C2025Q4R6 |\n| DUCC | Signed by every workspace member (Workbench requirement) |\n\n"
        "The IRB bundle stays local (STUDY-LIFECYCLE stage 5). No participant data in this packet.\n")


def main(all_tables, out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    figdir = out / "08_figure_pack"; figdir.mkdir(exist_ok=True)
    tabs = T.read_all(all_tables)
    import re
    hashes = sorted(set(re.findall(r"code ([0-9a-f]{12})", Path(all_tables).read_text())))
    code_hash = ", ".join(hashes)
    captions = []
    variable_definitions(out, code_hash)
    demographics(out, tabs, code_hash)
    for qn in ("Q1", "Q2"):
        if f"F04.flow {qn}" in tabs:
            flow_figure(tabs[f"F04.flow {qn}"], qn, out / f"03_cohort_flow_{qn}")
            for ext in ("png", "pdf"):
                (figdir / f"fig{1 if qn == 'Q1' else 2}_cohort_flow_{qn}.{ext}").write_bytes((out / f"03_cohort_flow_{qn}.{ext}").read_bytes())
            captions.append((f"fig{1 if qn == 'Q1' else 2}_cohort_flow_{qn}",
                             f"{qn} cohort construction: persons remaining after each protocol step and persons excluded. Masked: "
                             "steps of 1-20 persons are hidden and counted in the next shown step."))
    analysis_diagram(out / "04_analysis_diagram")
    for ext in ("png", "pdf"):
        (figdir / f"fig3_analysis_map.{ext}").write_bytes((out / f"04_analysis_diagram.{ext}").read_bytes())
    captions.append(("fig3_analysis_map", "Analysis map: the study questions, each prespecified analysis and its sensitivity "
                     "analyses (ANALYSIS-MAP), with the pipeline cell that runs it."))
    have = results(out, tabs, figdir, captions, code_hash)
    notebook(out)
    irb(out)
    (figdir / "CAPTIONS.md").write_text("# Figure pack: names and captions\n\n" + "\n".join(f"- **{a}**: {b}" for a, b in captions) + "\n")
    print("packet written to", out, "| results:", "yes" if have else "not yet (no F10 tables)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
