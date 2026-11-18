"""Merge the idea tables of Analysis/research/*.md into one ranking on a common score scale.

Usage: python tools/ideas_rank.py [--md OUT.md] [--csv OUT.csv]

Each research file has ideas tables whose rows start with `| <PREFIX>-A..` or `| <PREFIX>-B..`.
Task files (invsqrt, dqi, chemistry) give the gain as a % of the task cost; general.md gives it as a %
of the score. The global priority puts all ideas on the score scale:

    score_gain% = 100 * ((1 / (1 - task_gain/100)) ** (1/3) - 1)   (task files; about task_gain / 3)
    global_priority = P * score_gain% / effort_hours

Gain parsing takes the first number for "x (a-b)", the midpoint for ranges "a-b" / "a to b", and the
number after "≈" when present (general.md). Rows whose gain cannot be parsed keep the file's own
priority and are flagged. Cost: O(total lines).
"""
from __future__ import annotations

import argparse
import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = {"invsqrt": "INV", "dqi": "DQI", "chemistry": "CHEM", "general": "GEN"}
NUM = r"[-−+]?\d+(?:\.\d+)?"


def to_float(s: str) -> float:
    return float(s.replace("−", "-").replace("+", ""))


def parse_gain(cell: str) -> float | None:
    c = cell.replace("**", "").replace(",", "")
    m = re.search(r"≈\s*(" + NUM + r")\s*%", c)
    if m:
        return to_float(m.group(1))
    m = re.match(r"\s*(" + NUM + r")\s*\(", c)  # "0.67 (0.6-0.7)"
    if m:
        return to_float(m.group(1))
    m = re.match(r"\s*(" + NUM + r")\s*(?:–|-|to)\s*(" + NUM + r")", c)  # ranges
    if m:
        return (to_float(m.group(1)) + to_float(m.group(2))) / 2
    m = re.match(r"\s*(" + NUM + r")", c)
    if m:
        return to_float(m.group(1))
    m = re.search(r"(" + NUM + r")\s*%", c)  # e.g. "one invsqrt qubit = 0.77%"
    if m:
        return to_float(m.group(1))
    return None


def split_row(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def load() -> list[dict]:
    rows = []
    for name, prefix in FILES.items():
        path = os.path.join(ROOT, "Analysis", "research", f"{name}.md")
        if not os.path.exists(path):
            continue
        header = None
        for line in open(path, encoding="utf-8"):
            if line.startswith("| ID "):
                header = [h.lower() for h in split_row(line)]
                continue
            m = re.match(r"\|\s*(" + prefix + r"-[AB]\d+)\s*\|", line)
            if not (m and header):
                continue
            cells = split_row(line)
            rec = dict(zip(header, cells))
            gain_key = next(k for k in header if k.startswith("est"))
            p = to_float(rec["p"]) if "p" in rec else to_float(rec.get("p(success)", "nan"))
            effort_cell = next(v for k, v in rec.items() if k.startswith("effort"))
            effort = to_float(re.match(r"\s*(" + NUM + ")", effort_cell).group(1))
            g = parse_gain(rec[gain_key])
            own = rec.get("priority", "").replace("*", "")
            if g is None:
                score_gain, flag = None, "gain?"
            elif name == "general":
                score_gain, flag = g, ""
            else:
                score_gain = 100 * ((1 / (1 - g / 100)) ** (1 / 3) - 1) if g < 100 else g / 3
                flag = ""
            prio = p * score_gain / effort if score_gain is not None else None
            rows.append({
                "id": m.group(1), "task": name, "category": m.group(1).split("-")[1][0],
                "idea": rec.get("idea", ""), "gain_cell": rec[gain_key], "score_gain": score_gain,
                "p": p, "effort": effort, "priority": prio, "own_priority": own, "flag": flag,
            })
    return rows




EXCLUDED = {  # idea id -> reason (decisions recorded in Analysis/PLAN.md)
    "INV-A01": "excluded by decision S1.3 (track B only: no compiled third-party files); value carried by INV-A02..A04",
    "CHEM-A01": "excluded by decision S1.3 (track B only); value carried by CHEM-A03",
}
TASK_TITLES = {"invsqrt": "invsqrt", "dqi": "dqi", "chemistry": "Chemistry (dfthcsos)", "general": "Cross-cutting methods and tooling"}
CAT_TITLES = {"A": "A. Improve our code / tooling", "B": "B. Rewrite with a different technique / new pipeline"}


def md_cell(s: str) -> str:
    return s.replace("|", "/").replace("\n", " ")


def render_markdown(rows: list[dict], top: int = 30) -> str:
    """Global top list, then one table per task and category, all sorted by global priority."""
    def fmt(r, with_task):
        sg = f"{r['score_gain']:.2f}" if r["score_gain"] is not None else "?"
        pr = f"{r['priority']:.4f}" if r["priority"] is not None else "?"
        status = "EXCLUDED" if r["id"] in EXCLUDED else "TODO"
        cells = [r["id"]] + ([r["task"]] if with_task else []) + [md_cell(r["idea"]), md_cell(r["gain_cell"]), sg,
                 f"{r['p']:.2f}", f"{r['effort']:.0f}", pr, status]
        return "| " + " | ".join(cells) + " |"

    out = []
    live = [r for r in rows if r["id"] not in EXCLUDED]
    out.append(f"## 3. Global ranking (top {top}, all tasks)\n")
    out.append("| ID | Task | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |")
    out.append("|---|---|---|---|---:|---:|---:|---:|---|")
    out += [fmt(r, True) for r in live[:top]]
    out.append("")
    for i, (task, title) in enumerate(TASK_TITLES.items(), start=4):
        out.append(f"## {i}. {title}\n")
        out.append(f"Details (mechanism, risks, sources): `Analysis/research/{task}.md`.\n")
        for cat, ctitle in CAT_TITLES.items():
            sub = [r for r in rows if r["task"] == task and r["category"] == cat]
            if not sub:
                continue
            out.append(f"### {ctitle} ({len(sub)} ideas)\n")
            out.append("| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |")
            out.append("|---|---|---|---:|---:|---:|---:|---|")
            out += [fmt(r, False) for r in sub]
            out.append("")
    return "\n".join(out)


def write_markdown(path: str, top: int = 30) -> None:
    rows = load()
    rows.sort(key=lambda r: -(r["priority"] if r["priority"] is not None else -1))
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(render_markdown(rows, top))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv")
    ap.add_argument("--md", help="write the generated ranking tables (markdown) to this file")
    ap.add_argument("--top", type=int, default=30)
    a = ap.parse_args()
    rows = load()
    rows.sort(key=lambda r: -(r["priority"] if r["priority"] is not None else -1))
    for r in rows:
        sg = f"{r['score_gain']:.3f}" if r["score_gain"] is not None else "?"
        pr = f"{r['priority']:.4f}" if r["priority"] is not None else "?"
        print(f"{r['id']:10s} {r['task']:9s} sg={sg:>7s} P={r['p']:.2f} h={r['effort']:>5.0f} prio={pr:>8s} own={r['own_priority']:>7s} {r['flag']} | {r['idea'][:70]}")
    print(len(rows), "ideas")
    if a.md:
        write_markdown(a.md, a.top)
    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
