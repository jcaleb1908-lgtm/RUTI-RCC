"""Read the masked ALL_TABLES.txt into DataFrames (no Workbench access needed).

Each table in ALL_TABLES.txt is printed as
    --- <cell>.<id> <title> ---
    <pandas to_string(index=False) block, right-aligned>
    Note: ...
Columns are right-aligned, so each value ends where its header ends. Values are kept as text, because masked
cells read '<=20' or are blank.
"""
import re
from pathlib import Path

import pandas as pd

HEAD = re.compile(r"^--- (F\d\d\w*)\.(\S+) (.*) ---$")


def read_all(path):
    """Returns {"F04.flow Q1": {"title", "df", "note"}, ...}. A repeated id keeps the last copy (resumed runs)."""
    lines = Path(path).read_text().splitlines()
    out, i = {}, 0
    while i < len(lines):
        m = HEAD.match(lines[i])
        if not m:
            i += 1
            continue
        cell, tid, title = m.groups()
        body, j = [], i + 1
        while j < len(lines) and not lines[j].startswith("Note:") and not HEAD.match(lines[j]):
            if lines[j].strip():
                body.append(lines[j])
            j += 1
        note = lines[j][5:].strip() if j < len(lines) and lines[j].startswith("Note:") else ""
        df = parse_block(body)
        key = f"{cell}.{tid}"
        qm = re.match(r"(Q\d\w*)\b", title)
        if qm and tid in ("flow", "table1", "arms", "timing", "events_total"):
            key += " " + qm.group(1)
        out[key] = {"title": title, "df": df, "note": note}
        i = j + 1
    return out


def parse_block(body):
    if not body:
        return pd.DataFrame()
    head = body[0]
    spans = [(m.start(), m.end()) for m in re.finditer(r"\S+", head)]
    names = [head[a:b] for a, b in spans]
    if len(body) == 1 or len(names) == 0:
        return pd.DataFrame(columns=names)
    rows = []
    for ln in body[1:]:
        vals, prev = [], 0
        for k, (a, b) in enumerate(spans):
            if k == len(spans) - 1:
                vals.append(ln[prev:].strip())
            else:
                vals.append(ln[prev:b].strip())
                prev = b
        rows.append(vals)
    return pd.DataFrame(rows, columns=names)


def num(s):
    """Numeric where possible; masked or blank cells become NaN."""
    return pd.to_numeric(pd.Series(s).astype(str).str.replace(",", ""), errors="coerce")
