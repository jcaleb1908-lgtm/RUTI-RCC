"""RUN_IN_WORKBENCH.py - paste into one notebook cell in a JupyterLab app of the RUTI-RCC workspace.

Before the first run (WORKBENCH-2.0-GUIDE.md):
  1. The RUTI-RCC GitHub repo is a Git resource of the workspace, so it sits in ~/repos/<id>; `git pull` there.
  2. Add the CDR v9 data collection (C2025Q4R6). In a terminal: `wb resource list`, then `wb resolve --id=<cdr id>`.
  3. Create a workspace bucket for runs and raise the app autostop before a long run.
Never paste outputs back except runs/<run>/tables/ALL_TABLES.txt (masked aggregates). No AI tool inside the workspace.
"""
import hashlib, sys
from pathlib import Path

REPO = Path.home() / "repos" / "RUTI-RCC"          # <- adjust if the Git resource id differs
CDR = ""                                            # <- paste `wb resolve --id=<cdr id>` output, e.g. project.C2025Q4R6
BUCKET = ""                                         # <- workspace bucket resource id for the run folder
CELLS = ["F00_setup", "F01_inventory", "F02_extract", "F03_phenotypes", "F04_cohorts"]
RUN_CONTROL = {"new_run": True, "resume_from": None, "run_only": None}   # e.g. resume_from="F03" after a restart

P = REPO / "pipeline"
sys.path.insert(0, str(P))
# Forget modules from an earlier run in this kernel, so a `git pull` takes effect without a kernel restart.
for _m in ("aou_core", "_config", "rutirc", "cohort"):
    sys.modules.pop(_m, None)
import aou_core as C  # noqa: E402
CODE_HASH = hashlib.sha256(b"".join((P / f).read_bytes() for f in sorted(
    x.name for x in P.glob("*.py") if not x.name.startswith("test_")))).hexdigest()[:12]
C.CFG["cdr"] = CDR or (C.resolve_resource("cdr") or "")
if BUCKET:
    mnt = C.bucket_mount(BUCKET)
    if mnt is None:
        raise SystemExit(f"STOPPED: bucket {BUCKET} is not mounted (run `wb resource mount`).")
    C.CFG["run_root"] = str(mnt / "runs")
todo = CELLS
if RUN_CONTROL["resume_from"]:
    todo = CELLS[CELLS.index(RUN_CONTROL["resume_from"]):]
    RUN_CONTROL["new_run"] = False
if RUN_CONTROL["run_only"]:
    todo = [c for c in CELLS if c in RUN_CONTROL["run_only"]]
    RUN_CONTROL["new_run"] = "F00_setup" in todo
for cell in todo:
    print(f"\n==================== {cell} | code {CODE_HASH} ====================")
    C.CFG["footnote_extra"] = f"code {CODE_HASH}"
    exec(compile((P / f"{cell}.py").read_text(), cell, "exec"), {"RUN_CONTROL": RUN_CONTROL, "CODE_HASH": CODE_HASH,
                                                                   "__name__": "__main__"})
print("\nDone. Copy out ONLY:", C._run() / "tables" / "ALL_TABLES.txt")
