"""T10 seat 6: THE LANDING CHECK OF A TABLE STEP (the lead's 08:44Z order (3)), run on the box in the step's folder.

A table is fit to land in fin/ only if (1) table.json says the rows chain completed (table7c's round-2 refit found a
feasible point), (2) its final error_1 is at or under LIMIT1 (0.0200 Ha = the accuracy limit), and (3) the final table's
evaluation (vall_c.json, after the all-gates refit, C2's carrier gauge and the repair) has E_gap at or under GAPMAX (3.92;
the floor is 3.914925) and passes the accuracy test. Prints one line: CHECK_OK or CHECK_FAIL <reason> with the values, and
RETRY (a rows / refit failure: new tries may cure it) or NORETRY (the two-body part has no feasible weights).

usage: python3 chem_t10_tabcheck.py STEP_DIR [STEP_LOG]
"""
import json
import os
import sys

LIMIT1, GAPMAX = 0.0200, 3.92


def main():
    d = sys.argv[1]
    log = sys.argv[2] if len(sys.argv) > 2 else d.rstrip("/") + ".log"
    txt = open(log, errors="replace").read() if os.path.exists(log) else ""
    if "NO FEASIBLE WEIGHTS" in txt:
        print("CHECK_FAIL NORETRY: the free refit found no feasible weights (error_2 of the two-body part over the limit)")
        return
    try:
        t = json.load(open(os.path.join(d, "table.json")))
    except Exception as e:  # noqa: BLE001
        print(f"CHECK_FAIL RETRY: no table.json ({type(e).__name__})")
        return
    fin = t.get("final", {})
    e1, comp = float(fin.get("error_1", 9.9)), bool(t.get("complete"))
    try:
        v = json.load(open(os.path.join(d, "vall_c.json")))
        gap, lam, pa = float(v["E_gap"]), float(v["lambda_eff"]), bool(v.get("passes_accuracy"))
    except Exception:  # noqa: BLE001
        gap, lam, pa = 9.9, 99.0, False
    vals = (f"rows chain complete {comp}, table.json error_1 {1e3 * e1:.3f} mHa; final E_gap {gap:.6f}, "
            f"lambda_eff {lam:.6f}, passes_accuracy {pa}")
    why = []
    if not comp:
        why.append("the rows chain did not complete")
    if e1 > LIMIT1:
        why.append(f"error_1 over {1e3 * LIMIT1:.1f} mHa")
    if gap > GAPMAX:
        why.append(f"E_gap over {GAPMAX}")
    if not pa:
        why.append("fails the accuracy test")
    print(("CHECK_OK " if not why else "CHECK_FAIL RETRY: " + "; ".join(why) + ": ") + vals + price_note(d, lam))


def price_note(d, lam):
    """seat 6 (the lead's 10:27Z order): the table's price at ITS OWN shape by the hosts' code (tools/chem_t10_price.py),
    not the step's fixed 9 x 127 constants; quietly empty when no venv or table is at hand."""
    import subprocess
    tab = os.path.join(d, "all_c.npz")
    here = os.path.dirname(os.path.abspath(__file__))
    for py in ("~/qcc/.venv_work/bin/python", "~/qcc/.venv/bin/python"):
        py = os.path.expanduser(py)
        if os.path.exists(py) and os.path.exists(tab) and lam < 90:
            try:
                out = subprocess.run([py, "-I", os.path.join(here, "chem_t10_price.py"), tab, str(lam)],
                                     capture_output=True, text=True, timeout=120).stdout.strip()
                if out.startswith("PRICE"):
                    return "; " + out
            except Exception:  # noqa: BLE001
                pass
    return ""


if __name__ == "__main__":
    main()
