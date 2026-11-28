"""aou_core.py - reusable core for All of Us Researcher Workbench analyses (Controlled Tier, OMOP on BigQuery).

Generalized from the GDM-CVD V9 pipeline common block (jcaleb1908-lgtm/FINAL-GDM, 1_V9_Final_Analysis/pipeline/
_common_block.py, run 2026-10-07). Project-specific definitions (concept sets, windows, specs) do NOT live here:
they go in the project's own config block (see ANALYSIS-STANDARD.md, section 3).

What it gives every project
  * run folder + saved state (parquet) so each cell reloads its inputs after a kernel restart;
  * one BigQuery wrapper: import shim for the new Workbench (google-api-core too old), preflight, dry-run size,
    query parameters, the SQL text saved per query, and an offline override (DuckDB synthetic CDR in tests);
  * concept sets = roots + concept_ancestor descendants, each root name-checked against a regex;
  * disclosure: counts 1..HIDE_MAX shown as '<=20', complementary cells hidden (row partitions, groups,
    cross-tabs), masked cohort flow, cross-table leak check, person_id guard, percentiles instead of min/max;
  * cell skeleton: begin() / part() / end() - one failing part prints PART FAILED and the cell continues;
  * time helpers: completed age, first record inside a per-person window, record-in-any-window flag;
  * model helpers: Cox (statsmodels PHReg, delayed entry, robust SE), logistic, E-values, Rubin's rules.

Use in the Workbench (upload this file next to the notebook):
    import aou_core as C
    C.CFG.update(project="MYPROJ", version="MYPROJ-2026.10.08")
    C.begin("F00", "setup", create=True)
    ok = C.bq_preflight()
Offline (tests): C.CFG["q_override"] = fn(sql, params, label) -> DataFrame; C.CFG["run_dir"] = "<tmp folder>".
"""
import os
import re
import time
import uuid
import traceback
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIG (only generic settings; project definitions live in the project's config block)
# =============================================================================
CFG = {
    "project": "AOU",
    "version": "AOU-UNSET",                       # bump on every definition change; printed in every footnote
    "cdr": os.environ.get("WORKSPACE_CDR", ""),   # legacy variable; in Workbench 2.0 set it (project.dataset) or use resolve_resource()
    "seed": 20260101,
    "hide_max": 20,                               # All of Us DSDP: counts 1..20 never shown (confirm current policy)
    "run_root": "runs",                           # <cwd>/runs/<UTC stamp>_<hex>/
    "run_dir": None,                              # override (tests)
    "q_override": None,                           # override (tests): fn(sql, params, label) -> DataFrame
    "footnote_extra": "",                         # e.g. one-line spec text
}
STATE = {"cell": "", "run": None, "t0": time.time(), "tables": [], "failed": []}
_BQ = {}


def hide_max():
    return int(CFG["hide_max"])


# =============================================================================
# RUN FOLDER AND SAVED STATE
# =============================================================================
def run_dir(create=False):
    """Current run folder. create=True starts a new one and records it in LATEST.txt."""
    if CFG.get("run_dir"):
        p = Path(CFG["run_dir"]); p.mkdir(parents=True, exist_ok=True); return p
    base = Path.cwd() / CFG["run_root"]
    if create:
        p = base / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:6])
        p.mkdir(parents=True, exist_ok=False)
        (base / "LATEST.txt").write_text(p.name)
        return p
    if (base / "LATEST.txt").exists():
        p = base / (base / "LATEST.txt").read_text().strip()
        if p.exists():
            return p
    cand = sorted(d for d in base.glob("*") if d.is_dir()) if base.exists() else []
    if not cand:
        raise FileNotFoundError(f"{base} has no run folder: run the first cell with begin(..., create=True).")
    return cand[-1]


def _run():
    if STATE["run"] is None:
        STATE["run"] = run_dir()
    return STATE["run"]


def _state_path(name, internal=True):
    # Participant-level files carry _INTERNAL in the name: .gitignore and the copy-out rule key on it.
    return _run() / "state" / (f"{name}_INTERNAL.parquet" if internal else f"{name}.parquet")


def save_state(df, name, internal=True):
    p = _state_path(name, internal); p.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(p, index=False)
    except Exception:
        p = p.with_suffix(".pkl"); df.to_pickle(p)
    return p


def load_state(name, internal=True, required=True):
    p = _state_path(name, internal)
    if p.exists():
        return pd.read_parquet(p)
    if p.with_suffix(".pkl").exists():
        return pd.read_pickle(p.with_suffix(".pkl"))
    if required:
        raise FileNotFoundError(f"{p.name} missing in {_run()}: run the earlier cell that writes it.")
    return None


def has_state(name, internal=True):
    p = _state_path(name, internal)
    return p.exists() or p.with_suffix(".pkl").exists()


# =============================================================================
# BIGQUERY
# =============================================================================
class Prm:
    """Query parameter usable online (BigQuery) and offline (tests)."""
    def __init__(self, name, typ, value, array=False):
        self.name, self.typ, self.value, self.array = name, typ, value, array


def A(name, typ, values):
    """Array parameter. numpy scalars become Python values (the BigQuery client cannot serialize numpy int64)."""
    return Prm(name, typ, [v.item() if hasattr(v, "item") else v for v in values], True)


def S(name, typ, value):
    return Prm(name, typ, value.item() if hasattr(value, "item") else value, False)


BQ_FIX = ('In a new cell run:  %pip install --upgrade "google-api-core>=2.28"  then Kernel -> Restart and rerun.')


def _bq():
    """Import google.cloud.bigquery. In the new Workbench (2026-10) google-cloud-bigquery called helpers that the
    installed google-api-core did not have (AttributeError check_python_version): add them as no-ops and retry."""
    import sys
    try:
        from google.cloud import bigquery
        return bigquery
    except AttributeError as exc:
        if "check_python_version" not in str(exc) and "check_dependency_versions" not in str(exc):
            raise
    import google.api_core as _ac
    for fn in ("check_python_version", "check_dependency_versions"):
        if not hasattr(_ac, fn):
            setattr(_ac, fn, lambda *a, **k: None)
    for m in [m for m in list(sys.modules) if m.startswith("google.cloud.bigquery")]:
        del sys.modules[m]
    from google.cloud import bigquery
    print("note: google-api-core older than google-cloud-bigquery expects; shim applied. Permanent fix: " + BQ_FIX)
    return bigquery


def _client():
    if "client" not in _BQ:
        bigquery = _bq()
        proj = (os.environ.get("GOOGLE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT") or "").strip()
        _BQ["client"] = bigquery.Client(project=proj) if proj else bigquery.Client()
    return _BQ["client"]


def bq_preflight():
    """Import BigQuery, run SELECT 1 and check the CDR is set. Returns True when a full run can start."""
    if CFG.get("q_override"):
        return True
    try:
        list(_client().query("SELECT 1 AS ok").result(timeout=300))
    except Exception as exc:
        print(f"BigQuery preflight FAILED: {type(exc).__name__}: {str(exc)[:400]}\n{BQ_FIX}")
        return False
    if not CFG["cdr"]:
        print("BigQuery works but CFG['cdr'] is empty: set it to the CDR dataset (project.dataset) before any query.")
        return False
    print(f"BigQuery preflight: OK | CDR {CFG['cdr']} | project {os.environ.get('GOOGLE_CLOUD_PROJECT', '(unset)')}")
    rd = str(Path(CFG.get("run_dir") or (Path.cwd() / CFG["run_root"])).resolve())
    if "/workspace/" not in rd:
        print(f"note: run folder {rd} is on the app disk (lost if the app is deleted); a mounted workspace bucket "
              "(aou_core.bucket_mount) keeps it.")
    return True


def q(sql, prms=(), label="query", max_gib=None):
    """Run one query. Writes the SQL to <run>/sql/<label>.sql, prints the dry-run size, raises on failure.
    `{CDR}` in the SQL is replaced by CFG['cdr']. max_gib: refuse queries that would scan more."""
    t0 = time.time()
    sql = sql.replace("{CDR}", CFG["cdr"])
    try:
        d = _run() / "sql"; d.mkdir(parents=True, exist_ok=True)
        (d / (re.sub(r"[^A-Za-z0-9_]+", "_", label)[:60] + ".sql")).write_text(sql)
    except Exception:
        pass
    if CFG.get("q_override"):
        out = CFG["q_override"](sql, list(prms), label)
    else:
        bigquery = _bq()
        cl = _client()
        bp = [bigquery.ArrayQueryParameter(p.name, p.typ, p.value) if p.array else
              bigquery.ScalarQueryParameter(p.name, p.typ, p.value) for p in prms]
        dry = cl.query(sql, job_config=bigquery.QueryJobConfig(query_parameters=bp, dry_run=True, use_query_cache=False))
        gib = int(dry.total_bytes_processed or 0) / 1024 ** 3
        print(f"{label}: scans {gib:.1f} GiB ...", flush=True)
        if max_gib is not None and gib > max_gib:
            raise RuntimeError(f"{label} would scan {gib:.1f} GiB > max_gib {max_gib}")
        job = cl.query(sql, job_config=bigquery.QueryJobConfig(query_parameters=bp)).result(timeout=7200)
        try:
            out = job.to_dataframe(create_bqstorage_client=True)
        except Exception:
            out = job.to_dataframe(create_bqstorage_client=False)
    print(f"   {label}: {len(out):,} rows in {time.time() - t0:.0f} s", flush=True)
    return out


def resolve_resource(name):
    """Workbench 2.0: location of a workspace resource (bucket, BigQuery dataset, data collection reference).
    Tries the env variable WORKBENCH_<name> (not always pre-set), then `wb resource resolve --id=<name>`.
    Returns the string, or None when neither works (then set CFG['cdr'] by hand from the workspace Resources tab)."""
    import subprocess
    v = os.environ.get(f"WORKBENCH_{name}", "").strip()
    if v:
        return v
    try:
        r = subprocess.run(["wb", "resource", "resolve", f"--id={name}"], capture_output=True, text=True, timeout=60)
        out = r.stdout.strip()
        return out or None
    except Exception:
        return None


def bucket_mount(name, home=None):
    """Workbench 2.0 mounts workspace buckets under $HOME/workspace/<resource folders>/<resource id> (Cloud Storage FUSE).
    Returns the mounted directory of bucket resource `name`, or None. Put the run folder there
    (CFG['run_root'] = str(bucket_mount('my_bucket') / 'runs')): files on the app disk are lost when the app is deleted."""
    base = Path(home or os.path.expanduser("~")) / "workspace"
    if not base.exists():
        return None
    # Look only at base/<name> and base/<folder>/<name>. Never rglob: ~/workspace also mounts the data buckets
    # (Cloud Storage FUSE), and walking them can take hours.
    for cand in [base / name] + sorted(d / name for d in base.iterdir() if d.is_dir()):
        if cand.is_dir():
            return cand
    return None


def gcs_cp(src, dst, billing_project=None, timeout=1800):
    """Copy with `gcloud storage cp` (gsutil is deprecated and leaves the Google Cloud CLI bundle in March 2027).
    billing_project: for requester-pays buckets (defaults to GOOGLE_CLOUD_PROJECT / GOOGLE_PROJECT)."""
    import subprocess
    proj = billing_project or os.environ.get("GOOGLE_CLOUD_PROJECT") or os.environ.get("GOOGLE_PROJECT")
    cmd = ["gcloud", "storage", "cp", src, dst] + ([f"--billing-project={proj}"] if proj else [])
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"gcloud storage cp failed ({r.returncode}): {r.stderr.strip()[:400]}")
    return dst


def to_date(s):
    """Any date-like column (incl. BigQuery dbdate, tz-aware) -> naive datetime64[ns]."""
    try:
        s = pd.to_datetime(pd.Series(s), errors="coerce", utc=True)
    except (TypeError, ValueError):
        s = pd.to_datetime(pd.Series(s).astype(str), errors="coerce", utc=True)
    return s.dt.tz_localize(None).astype("datetime64[ns]")


# =============================================================================
# CONCEPT SETS
# =============================================================================
ROOT_SQL = """SELECT concept_id, concept_name, domain_id, vocabulary_id, concept_class_id, standard_concept, invalid_reason
FROM `{CDR}.concept` WHERE concept_id IN UNNEST(@ids)"""
DESC_SQL = """SELECT ca.ancestor_concept_id AS root_id, c.concept_id, c.concept_name, c.domain_id, c.vocabulary_id,
       c.standard_concept, ca.min_levels_of_separation AS levels
FROM `{CDR}.concept_ancestor` ca JOIN `{CDR}.concept` c ON c.concept_id = ca.descendant_concept_id
WHERE ca.ancestor_concept_id IN UNNEST(@ids)"""
SRC_MAP_SQL = """SELECT s.concept_id AS source_concept_id, s.vocabulary_id AS source_vocabulary, s.concept_code AS source_code,
       s.concept_name AS source_name, r.concept_id_2 AS standard_concept_id
FROM `{CDR}.concept` s
LEFT JOIN `{CDR}.concept_relationship` r ON r.concept_id_1 = s.concept_id AND r.relationship_id = 'Maps to'
WHERE s.vocabulary_id IN ('ICD9CM', 'ICD10CM') AND REGEXP_CONTAINS(s.concept_code, @rx)"""


def expand_concept_set(key, roots, intent_rx=None, keep_failed=False):
    """roots: [(concept_id, expected-name regex)]. Returns (concepts, root_check).
    concepts: every descendant (incl. the root itself) with set_key, root_id, levels, outside_intent.
    A root whose name does not match its regex is reported (name_ok = False) and dropped unless keep_failed.
    A non-standard root has no descendants in concept_ancestor: it is reported (n_descendants = 0)."""
    ids = [int(r[0]) for r in roots]
    R = q(ROOT_SQL, [A("ids", "INT64", ids)], f"{key} roots")
    rx = {int(r[0]): r[1] for r in roots}
    R["expected_rx"] = R["concept_id"].astype(int).map(rx)
    R["name_ok"] = [bool(re.search(e, str(n), re.I)) for e, n in zip(R["expected_rx"], R["concept_name"])]
    missing = sorted(set(ids) - set(R["concept_id"].astype(int)))
    if missing:
        R = pd.concat([R, pd.DataFrame({"concept_id": missing, "concept_name": "NOT IN CDR", "name_ok": False})],
                      ignore_index=True)
    R["name_ok"] = R["name_ok"].astype(bool)
    use = [int(c) for c, ok in zip(R["concept_id"], R["name_ok"]) if ok or keep_failed]
    for _, r in R.loc[~R["name_ok"]].iterrows():
        print(f"ROOT NAME CHECK FAILED {key}: {r['concept_id']} '{r['concept_name']}' !~ /{r.get('expected_rx')}/"
              + (" (kept)" if keep_failed else " (dropped)"))
    D = q(DESC_SQL, [A("ids", "INT64", use)], f"{key} descendants") if use else pd.DataFrame(
        columns=["root_id", "concept_id", "concept_name", "domain_id", "vocabulary_id", "standard_concept", "levels"])
    D = D.sort_values("levels").drop_duplicates("concept_id")
    D.insert(0, "set_key", key)
    D["outside_intent"] = (~D["concept_name"].astype(str).str.contains(intent_rx, case=False, regex=True)
                           if intent_rx else False)
    nd = D.groupby("root_id").size()
    R["n_descendants"] = R["concept_id"].astype(int).map(nd).fillna(0).astype(int)
    R.insert(0, "set_key", key)
    return D.reset_index(drop=True), R


def source_code_map(rx):
    """ICD-9-CM / ICD-10-CM source codes matching rx and the standard concept each maps to. Use it to check that the
    source codes you mean (e.g. ICD-9 648.8x) land inside your concept set; codes that map outside it are missed."""
    return q(SRC_MAP_SQL, [S("rx", "STRING", rx)], f"source map {rx}")


# =============================================================================
# PROGRAM TABLES (All of Us Support Hub)
# =============================================================================
# Enrollment proxy: first observation under the "Consent PII" module (primary consent, not EHR consent).
PRIMARY_CONSENT_SQL = """SELECT o.person_id, MIN(o.observation_date) AS primary_consent_date
FROM `{CDR}.concept` c
JOIN `{CDR}.concept_ancestor` ca ON ca.ancestor_concept_id = c.concept_id
JOIN `{CDR}.observation` o ON o.observation_concept_id = ca.descendant_concept_id
WHERE c.concept_name = 'Consent PII' AND c.concept_class_id = 'Module'
GROUP BY o.person_id"""
# aou_death holds EHR and HealthPro rows, several per person possible; the primary flag keeps one.
PRIMARY_DEATH_SQL = """SELECT person_id, death_date
FROM `{CDR}.aou_death` WHERE primary_death_record = TRUE"""


def enrollment_dates():
    """person_id -> primary_consent_date (delayed entry for left truncation at enrollment, AL5)."""
    d = q(PRIMARY_CONSENT_SQL, [], "primary consent date")
    d["primary_consent_date"] = to_date(d["primary_consent_date"])
    return d.set_index("person_id")["primary_consent_date"]


def death_dates():
    """person_id -> death_date from the primary death record. Stops if the flag is not unique per person."""
    d = q(PRIMARY_DEATH_SQL, [], "primary death record")
    if d["person_id"].duplicated().any():
        raise ValueError("aou_death: more than one primary_death_record per person")
    d["death_date"] = to_date(d["death_date"])
    return d.set_index("person_id")["death_date"]


# =============================================================================
# DISCLOSURE (All of Us Data and Statistics Dissemination Policy)
# =============================================================================
COUNT_RX = r"(^n$|^n_|_n$|^events?$|^events_|_events$|^count$|^women|_women$|^persons?$|_persons$|^n[0-9])"
EST_RX = r"(^OR|^HR|^RR|^RD|^RERI|^estimate|^lo$|^hi$|^p$|^p_|95% CI|^pct|^rate|^share|^CIF|^E_value)"
DEFAULT_PARTITIONS = [("n", "n_exposed", "n_unexposed"), ("events", "events_exposed", "events_unexposed")]


def _small(v):
    try:
        x = float(v)
    except Exception:
        return False
    return bool(np.isfinite(x)) and 1 <= x <= hide_max()


def suppress(df, counts=None, partitions=None, pct=None, group_cols=None, cross=False, deps=None, group_all=False):
    """Copy of df with counts 1..HIDE_MAX -> '<=20' and complementary cells hidden.
    partitions: tuples (total, part1, part2, ...) on one row; if one part is hidden, or hidden parts add up to 1..20, or
      total minus visible parts is 1..20, the smallest other part is hidden too.
    group_cols: rows of a group partition a total (Table 1 levels); [] = whole table is one group. group_all=True hides
      the whole group column when any row is hidden (overlapping categories).
    cross: also apply the rule across the count columns of each row (cross-tabs).
    pct / deps: {derived column: count column(s)} blanked when a source count is hidden."""
    out = df.copy()
    if counts is None:
        counts = [c for c in out.columns if re.search(COUNT_RX, str(c))]
    counts = [c for c in counts if c in out.columns]
    H = pd.DataFrame(False, index=out.index, columns=counts)
    V = pd.DataFrame({c: pd.to_numeric(out[c], errors="coerce").astype(float) for c in counts}, index=out.index)
    for c in counts:
        H[c] = out[c].apply(_small)
    parts = [tuple(c for c in p if c in counts) for p in (partitions if partitions is not None else DEFAULT_PARTITIONS)]
    parts = [p for p in parts if len(p) >= 2]

    def hide_smallest(i, cols):
        cand = [c for c in cols if not H.at[i, c] and pd.notna(V.at[i, c]) and V.at[i, c] > 0]
        if cand:
            H.at[i, min(cand, key=lambda c: V.at[i, c])] = True

    def hidden_sum_small(h, v):
        return int(h.sum()) >= 1 and _small(v[h].fillna(0).sum())

    if cross and group_cols is None:
        group_cols = []
    if group_cols is not None:
        gkeys = (pd.Series("all", index=out.index) if len(group_cols) == 0 else
                 out[group_cols].astype(str).agg("|".join, axis=1))
    for _ in range(10):
        before = int(H.values.sum())
        for p in parts:
            tot, pp = p[0], list(p[1:])
            for i in out.index:
                k = int(H.loc[i, pp].sum())
                if H.at[i, tot] or (len(pp) < 2 and k):
                    continue
                if (k == 1 and len(pp) > 1) or (k >= 2 and hidden_sum_small(H.loc[i, pp], V.loc[i, pp])):
                    hide_smallest(i, pp)
                elif pd.notna(V.at[i, tot]):
                    rem = V.at[i, tot] - V.loc[i, pp][~H.loc[i, pp]].fillna(0).sum()
                    if _small(rem):
                        hide_smallest(i, pp)
        if cross:
            for i in out.index:
                k = int(H.loc[i].sum())
                if k == 1 or (k >= 2 and hidden_sum_small(H.loc[i], V.loc[i])):
                    hide_smallest(i, counts)
        if group_cols is not None:
            for c in counts:
                for _, idx in gkeys.groupby(gkeys).groups.items():
                    k = int(H.loc[idx, c].sum())
                    if group_all and k >= 1:
                        H.loc[idx, c] = H.loc[idx, c] | V.loc[idx, c].gt(0)
                    elif len(idx) >= 2 and (k == 1 or (k >= 2 and hidden_sum_small(H.loc[idx, c], V.loc[idx, c]))):
                        cand = V.loc[idx, c][~H.loc[idx, c] & V.loc[idx, c].gt(0)]
                        if len(cand):
                            H.at[cand.idxmin(), c] = True
        if int(H.values.sum()) == before:
            break
    for c in counts:
        isint = V[c].notna() & np.isclose(V[c].fillna(0) % 1, 0)
        out[c] = out[c].astype(object)
        out.loc[isint, c] = V.loc[isint, c].astype("int64")
        out.loc[H[c], c] = "<=20"
    dd = dict(deps or {})
    dd.update({k: [v] for k, v in (pct or {}).items()})
    for dc, srcs in dd.items():
        srcs = [s for s in ([srcs] if isinstance(srcs, str) else srcs) if s in H.columns]
        if dc in out.columns and srcs:
            m = H[srcs].any(axis=1)
            out[dc] = out[dc].astype(object)
            out.loc[m, dc] = "<=20" if dc in (pct or {}) else ""
    return out


def leak_hidden(M, vec, known=()):
    """Cross-row / cross-table check. M: {key: boolean mask over persons}; vec: boolean (the counted event or 1 for
    persons). A count 1..20 (a key, a membership atom, or the sum of hidden levels of one grouping - tuple keys share a
    grouping by their first element) is disclosed when it is a linear combination of shown keys. Shown keys (except
    `known`, already shown elsewhere) are hidden, smallest first, until nothing small is computable.
    Returns (hidden keys, counts)."""
    hm = hide_max()
    keys = list(M)
    A_ = np.column_stack([np.asarray(M[k], bool) for k in keys])
    pat, inv = np.unique(A_, axis=0, return_inverse=True)
    ac = np.bincount(inv.ravel(), weights=np.asarray(vec, float), minlength=len(pat))
    inc = pat.T.astype(float)
    cnt = {k: int(round(inc[i] @ ac)) for i, k in enumerate(keys)}
    hidden = {k for k in keys if 1 <= cnt[k] <= hm}
    small_atoms = np.where((ac >= 1) & (ac <= hm) & pat.any(axis=1))[0]
    hid_rows = [keys.index(k) for k in hidden]
    groups = {}
    for i, k in enumerate(keys):
        if isinstance(k, tuple) and len(k) > 1:
            groups.setdefault(k[0], []).append(i)

    def basis(vis):
        if not vis:
            return np.zeros((0, inc.shape[1]))
        _, sv, vt = np.linalg.svd(inc[[keys.index(k) for k in vis]], full_matrices=False)
        return vt[sv > 1e-9 * max(sv.max(), 1.0)]

    def in_span(v, Q):
        return len(Q) > 0 and np.linalg.norm(v - (Q @ v) @ Q) < 1e-6

    def first_disclosed(vis):
        Q = basis(vis)
        if not len(Q):
            return None
        if len(small_atoms):
            res = 1.0 - (Q[:, small_atoms] ** 2).sum(0)
            hit = np.where(res < 1e-9)[0]
            if len(hit):
                v = np.zeros(inc.shape[1]); v[small_atoms[hit[0]]] = 1.0
                return v
        for r in hid_rows:
            if in_span(inc[r], Q):
                return inc[r]
        for g in groups.values():
            h = [i for i in g if keys[i] in hidden]
            if len(h) >= 2:
                v = inc[h].sum(0)
                if 1 <= round(v @ ac) <= hm and in_span(v, Q):
                    return v
        return None

    for _ in range(len(keys)):
        vis = [k for k in keys if k not in hidden]
        v = first_disclosed(vis)
        if v is None:
            break
        cand = sorted([k for k in vis if cnt[k] > 0 and k not in known], key=lambda k: cnt[k])
        if not cand:
            break
        hidden.add(next((k for k in cand if not in_span(v, basis([x for x in vis if x != k]))), cand[0]))
    return hidden, cnt


def flow_mask(steps):
    """Cohort flow [(label, n), ...] -> DataFrame. A step removing 1..20 persons is not shown (n blank); its exclusion
    is carried into the next shown step, so no small exclusion is recoverable by subtraction."""
    out, last = [], None
    for label, n in steps:
        d = None if last is None else last - n
        if d is not None and 1 <= d <= hide_max():
            out.append({"step": label, "n": np.nan, "excluded": "<=20 (counted in the next shown step)"})
        else:
            out.append({"step": label, "n": n, "excluded": "" if d is None else d})
            last = n
    return pd.DataFrame(out)


def pctl(s, ps=(25, 50, 75)):
    """Percentiles as 'p25 / p50 / p75' (min/max are never shown: they are single persons)."""
    s = pd.to_numeric(pd.Series(s), errors="coerce").dropna()
    if len(s) == 0:
        return ""
    return " / ".join(f"{x:.1f}" for x in np.percentile(s, ps))


def footnote(note=None):
    parts = [CFG["footnote_extra"]] if CFG.get("footnote_extra") else []
    if note:
        parts.append(note)
    parts.append(f"Definitions version {CFG['version']}; CDR {CFG['cdr'] or 'unset'}; counts 1-{hide_max()} shown as <=20 "
                 "(and complementary cells hidden).")
    return "Note: " + " | ".join(parts)


def _fmt(df):
    o = df.copy()
    for c in o.columns:
        if o[c].dtype == object or pd.api.types.is_string_dtype(o[c]):
            o[c] = o[c].astype(object).where(o[c].notna(), "")
        elif pd.api.types.is_float_dtype(o[c]):
            o[c] = o[c].map(lambda v: "" if pd.isna(v) else (f"{v:.3f}" if abs(v) < 1000 else f"{v:,.0f}"))
    return o


def show(tid, title, df, counts=None, partitions=None, pct=None, group_cols=None, cross=False, deps=None,
         group_all=False, note=None, max_rows=60):
    """Suppress, save (CSV + _index.csv registry + <cell>_tables.txt) and print one aggregate table.
    Estimates in a row with 1..20 exposed are blanked. A person_id column stops the cell."""
    full = f"{STATE['cell']}.{tid} {title}"
    print(f"\n--- {full} ---")
    foot = footnote(note)
    if df is None or len(df) == 0:
        print("(no rows)"); print(foot); return None
    if any(re.search(r"person_id|research_id", str(c), re.I) for c in df.columns):
        raise RuntimeError("person-level identifier in an aggregate table: refused")
    sv = suppress(df.reset_index(drop=True), counts, partitions, pct, group_cols, cross, deps, group_all)
    if "n_exposed" in df.columns:
        few = df["n_exposed"].reset_index(drop=True).apply(_small)
        est = [c for c in sv.columns if re.search(EST_RX, str(c))]
        if few.any() and est:
            sv[est] = sv[est].astype(object)
            sv.loc[few, est] = ""
            sv["note"] = sv["note"].astype(object) if "note" in sv.columns else ""
            sv.loc[few, "note"] = "estimate not shown (<=20 exposed)"
    tdir = _run() / "tables"; tdir.mkdir(parents=True, exist_ok=True)
    fname = f"{STATE['cell']}_{tid}.csv"
    sv.to_csv(tdir / fname, index=False)
    reg = tdir / "_index.csv"
    row = pd.DataFrame([{"cell": STATE["cell"], "tid": str(tid), "title": title, "footnote": foot, "file": fname,
                         "rows": len(sv), "written_utc": datetime.now(timezone.utc).isoformat()}])
    if reg.exists():
        old = pd.read_csv(reg, dtype={"tid": str})
        old = old.loc[~(old["cell"].eq(STATE["cell"]) & old["tid"].eq(str(tid)))]
        row = pd.concat([old, row], ignore_index=True)
    row.to_csv(reg, index=False)
    txt = _fmt(sv).head(max_rows).to_string(index=False)
    if len(sv) > max_rows:
        txt += f"\n(showing {max_rows} of {len(sv)} rows; full table: tables/{fname})"
    print(txt); print(foot, flush=True)
    with open(tdir / "ALL_TABLES.txt", "a") as fh:
        fh.write(f"\n--- {full} ---\n{_fmt(sv).to_string(index=False)}\n{foot}\n")
    STATE["tables"].append(fname)
    return sv


# =============================================================================
# CELL SKELETON
# =============================================================================
def begin(cell, title, create=False):
    STATE.update(cell=cell, t0=time.time(), tables=[])
    STATE["run"] = run_dir(create=create)
    (STATE["run"] / "tables").mkdir(parents=True, exist_ok=True)
    print(f"##### BEGIN {cell} | {title} #####", flush=True)
    print(f"run folder: {STATE['run']} | definitions {CFG['version']} | CDR {CFG['cdr'] or '(unset)'}", flush=True)


def part(name, fn):
    """Run one part in try/except: a failure prints PART FAILED with the traceback lines from your code."""
    t0 = time.time()
    print(f"\n========== {STATE['cell']} | PART {name} ==========", flush=True)
    if fn.__doc__:
        print("Purpose: " + " ".join(fn.__doc__.split()), flush=True)
    try:
        fn()
    except Exception as exc:
        print(f"PART FAILED: {type(exc).__name__}: {str(exc)[:300]}")
        tb = [f"  line {fr.lineno} in {fr.name}: {fr.line}" for fr in traceback.extract_tb(exc.__traceback__)
              if "-packages" not in fr.filename]
        print("\n".join(tb[-6:]))
        STATE["failed"].append(f"{STATE['cell']} part {name}")
    print(f"[{name.split(' ')[0]} {time.time() - t0:.0f} s; cell total {time.time() - STATE['t0']:.0f} s]", flush=True)


def end():
    print(f"\ntables written by {STATE['cell']}: {len(STATE['tables'])}")
    if STATE["failed"]:
        print("FAILED PARTS so far: " + "; ".join(STATE["failed"]))
    print(f"##### END {STATE['cell']} #####", flush=True)


# =============================================================================
# TIME HELPERS
# =============================================================================
def completed_age(birth, at):
    """Age in completed years (not the calendar-year difference)."""
    birth, at = to_date(birth).reset_index(drop=True), to_date(at).reset_index(drop=True)
    before_bday = (at.dt.month * 100 + at.dt.day) < (birth.dt.month * 100 + birth.dt.day)
    return (at.dt.year - birth.dt.year - before_bday.astype(int)).astype("Int64")


def first_in_window(rec, F, lo=None, hi=None, lo_incl=True, hi_incl=True, date_col="date"):
    """Per-person first record date with lo <(=) date <(=) hi. rec: person_id, date. F: person_id + the window
    columns named by lo / hi (None = open). Returns a Series indexed by person_id (persons without a record absent)."""
    cols = ["person_id"] + [c for c in (lo, hi) if c]
    m = rec[["person_id", date_col]].merge(F[cols], on="person_id", how="inner")
    k = pd.Series(True, index=m.index)
    if lo:
        k &= (m[date_col] >= m[lo]) if lo_incl else (m[date_col] > m[lo])
    if hi:
        k &= (m[date_col] <= m[hi]) if hi_incl else (m[date_col] < m[hi])
    return m.loc[k].groupby("person_id")[date_col].min()


def in_windows(rec, W, date_col="date"):
    """Boolean per record: date inside any [ws, we] window of the same person (W: person_id, ws, we)."""
    if rec is None or len(rec) == 0 or W is None or len(W) == 0:
        return pd.Series(False, index=getattr(rec, "index", None))
    r = rec[["person_id", date_col]].reset_index().rename(columns={"index": "_rid"})
    m = r.merge(W[["person_id", "ws", "we"]], on="person_id", how="inner")
    hit = m.loc[(m[date_col] >= m["ws"]) & (m[date_col] <= m["we"]), "_rid"].unique()
    return pd.Series(rec.index.isin(hit), index=rec.index)


# =============================================================================
# MODELS
# =============================================================================
def design(d, covs, cat_ref=None, numeric=()):
    """Design matrix: numeric covariates as is, others one-hot with the reference level dropped
    (cat_ref[col], else the most frequent level). Returns (X, term names)."""
    cat_ref = cat_ref or {}
    cols = []
    for c in covs:
        if c in numeric or (pd.api.types.is_numeric_dtype(d[c]) and c not in cat_ref):
            cols.append(pd.to_numeric(d[c], errors="coerce").astype(float).rename(c))
        else:
            v = d[c].astype(str)
            ref = cat_ref.get(c, v.value_counts().index[0])
            for lev in sorted(set(v) - {ref}):
                cols.append((v == lev).astype(float).rename(f"{c}[{lev}]"))
    X = pd.concat(cols, axis=1) if cols else pd.DataFrame(index=d.index)
    return X, list(X.columns)


DIVERGENT_LOG = 10.0   # |log HR| or |log OR| above this on any term: the fit is not estimable (M12)


def _coef_sheet(names, params, bse, pvalues, ratio):
    """Every term of a fitted model: beta, SE, ratio (HR or OR), 95% CI, p, and a divergence flag."""
    b, se, p = (np.asarray(v, dtype=float) for v in (params, bse, pvalues))
    sheet = pd.DataFrame({"term": list(names), "beta": b, "se": se, ratio: np.exp(b),
                          "lo": np.exp(b - 1.96 * se), "hi": np.exp(b + 1.96 * se), "p": p})
    sheet["divergent"] = ~np.isfinite(b) | ~np.isfinite(se) | (np.abs(b) > DIVERGENT_LOG)
    return sheet


def _result(sheet, x, ratio, n, events, caught):
    """Exposure row (keys unchanged) + the full coefficient sheet + status. status is
    'not_estimable:divergent_coefficient' when ANY term diverged, even if the exposure row looks normal."""
    term = [t for t in sheet["term"] if t == x or t.startswith(f"{x}[")][0]
    row = sheet.loc[sheet["term"] == term].iloc[0]
    bad = sheet.loc[sheet["divergent"], "term"].tolist()
    return {"term": term, ratio: float(row[ratio]), "lo": float(row["lo"]), "hi": float(row["hi"]),
            "p": float(row["p"]), "n": int(n), "events": int(events), "beta": float(row["beta"]),
            "se": float(row["se"]),
            "status": "not_estimable:divergent_coefficient" if bad else "estimated",
            "divergent_terms": bad, "warnings": sorted({str(w.message) for w in caught}),
            "coefs": sheet}


def fit_cox(d, t, y, x, covs=(), entry=None, cat_ref=None, numeric=(), robust=True, ties="efron"):
    """Cox model (statsmodels PHReg). Returns the exposure term x (HR, lo, hi, p, n, events, beta, se),
    the full coefficient sheet ("coefs", every term), "status" and "divergent_terms". Read "status"
    before reporting: one diverged nuisance term makes the whole fit not estimable (M12).
    entry: delayed-entry column (left truncation, e.g. time from time zero to enrollment)."""
    import statsmodels.api as sm
    X, names = design(d, [x] + list(covs), cat_ref, numeric)
    ok = X.notna().all(axis=1) & d[t].notna() & d[y].notna()
    m = sm.PHReg(d.loc[ok, t].astype(float), X.loc[ok], status=d.loc[ok, y].astype(int),
                 entry=None if entry is None else d.loc[ok, entry].astype(float), ties=ties)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = m.fit(groups=np.arange(int(ok.sum())) if robust else None)   # groups -> robust (sandwich) SE
    sheet = _coef_sheet(names, r.params, r.bse, r.pvalues, "HR")
    return _result(sheet, x, "HR", ok.sum(), d.loc[ok, y].sum(), caught)


def fit_logit(d, y, x, covs=(), cat_ref=None, numeric=()):
    """Logistic model. Returns the exposure term (OR, lo, hi, p, n, events, beta, se), the full
    coefficient sheet without the intercept ("coefs"), "status" and "divergent_terms" (M12)."""
    import statsmodels.api as sm
    X, names = design(d, [x] + list(covs), cat_ref, numeric)
    ok = X.notna().all(axis=1) & d[y].notna()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = sm.GLM(d.loc[ok, y].astype(float), sm.add_constant(X.loc[ok], has_constant="add"),
                   family=sm.families.Binomial()).fit()
    sheet = _coef_sheet(names, r.params[names], r.bse[names], r.pvalues[names], "OR")
    return _result(sheet, x, "OR", ok.sum(), d.loc[ok, y].sum(), caught)


def evalue_rr(rr):
    """E-value for a risk ratio (VanderWeele & Ding 2017)."""
    if rr is None or not np.isfinite(rr) or rr <= 0:
        return np.nan
    rr = 1 / rr if rr < 1 else rr
    return rr + np.sqrt(rr * (rr - 1))


def hr_to_rr_common(hr):
    """HR -> approximate RR when the outcome is common (> 15% by end of follow-up)."""
    return (1 - 0.5 ** np.sqrt(hr)) / (1 - 0.5 ** np.sqrt(1 / hr))


def evalue_ci(est, lo, hi):
    """E-value for the CI limit closest to the null (1 if the CI includes 1)."""
    if lo <= 1 <= hi:
        return 1.0
    return evalue_rr(lo if est > 1 else hi)


def rubin(qs, us):
    """Rubin's rules: estimates qs and within-imputation variances us -> (estimate, se, df)."""
    qs, us = np.asarray(qs, float), np.asarray(us, float)
    m = len(qs)
    qbar, ubar, b = qs.mean(), us.mean(), qs.var(ddof=1)
    tv = ubar + (1 + 1 / m) * b
    r = (1 + 1 / m) * b / ubar if ubar > 0 else np.inf
    dfree = (m - 1) * (1 + 1 / r) ** 2 if np.isfinite(r) and r > 0 else np.inf
    return qbar, np.sqrt(tv), dfree
