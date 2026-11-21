"""chem_ld2_hb: the number of 16-address high blocks for the C-4 tools (builder C2, tenth holder).

HB = env CHEM_HB (default 77 = the 1,220-row tables: 76 + 8 x 143). A table at B rows per group has
ceil((76 + 8 B) / 16) blocks (B 135: 73, B 129: 70); the B-ladder scripts export CHEM_HB for it. With CHEM_HB unset
every C-4 tool behaves exactly as before (77 blocks, last block 76, padding rows 1220..1231).
"""
import os

HB = int(os.environ.get("CHEM_HB", "77"))
LAST = HB - 1              # index of the last (partial) block
NFULL = 16 * HB            # rows of the padded table


def check_rows(A):
    """A table of A rows must end inside the last block."""
    assert NFULL - 16 < A <= NFULL, f"table has {A} rows but CHEM_HB={HB} expects {NFULL - 15}..{NFULL}"
