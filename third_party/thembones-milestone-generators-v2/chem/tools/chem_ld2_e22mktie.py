"""C2 (22nd holder): writes tools/chem_ld2_e22tie.py = a pinned COPY of R59-CHEM-GROUPS' tools/r59_chem_groups_tie.py
(read-only there) with one change: the tie group's relation count c_tot = len(forced) + bon - 2 takes
C2_TIE_EXTRA more free relations (env, default 0 = R59's count). A probe: with a row kill on group 2 (dim W 74, 16
relations available) and 3 ties, EXTRA 2 asks for 16 relations = 37 products instead of 38.
  python3 tools/chem_ld2_e22mktie.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = "            c_tot = min(len(rel), len(fmask) + bon - 2)\n"
NEW = ("            c_tot = min(len(rel), len(fmask) + bon - 2 + int(os.environ.get(\"C2_TIE_EXTRA\", \"0\")))"
       "   # C2 copy: + EXTRA\n")


def main():
    src = open(os.path.join(HERE, "r59_chem_groups_tie.py")).read()
    assert src.count(OLD) == 1, "upstream changed: re-read r59_chem_groups_tie.py"
    out = src.replace(OLD, NEW)
    if "\nimport os\n" not in out:
        out = out.replace("\nimport numpy as np\n", "\nimport os\nimport numpy as np\n", 1)
    out = "# C2 (22nd holder): PINNED COPY of tools/r59_chem_groups_tie.py made by tools/chem_ld2_e22mktie.py (C2_TIE_EXTRA)\n" + out
    open(os.path.join(HERE, "chem_ld2_e22tie.py"), "w").write(out)
    print("written chem_ld2_e22tie.py; import os present:", "\nimport os\n" in out)


if __name__ == "__main__":
    main()
