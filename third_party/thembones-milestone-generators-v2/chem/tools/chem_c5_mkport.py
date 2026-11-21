#!/usr/bin/env python3
"""chem_c5_mkport.py: seat C5 (Nov 19 20:12Z). Writes tools/sh/chem_c5_port.sh = builder C4's tools/sh/chem_c4_port5.sh with ONE
change: the EXACT path (FAST=0, = port4's path) builds a WIRED k0 in the port by loader2 on the K0C cores, as port5's FAST
path does (port4 had no in-port k0, so auto6 called every k0-wired layout 'no plan'). Run from the tools folder.
    python3 chem_c5_mkport.py
"""
import hashlib
import os

here = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(here, "sh", "chem_c4_port5.sh")
src = open(src_path).read()
old = "else\n( CA=(${RA[@]}); NC=${#CA[@]}; right_chain ) &\nPR=$!\n"
assert src.count(old) == 1, src.count(old)
new = ("else\n"
       "# seat C5 (Nov 19 20:12Z, the lead's 20:00Z: the exact port of a k0-WIRED table): k0 by loader2 here too, as the FAST\n"
       "# path does (port4 had no in-port k0); batch 0 is toggled by nobody (k1 plain), checked at the end as in the FAST path\n"
       "if [ \"$(ngof 0)\" != 30 ]; then\n"
       "  K0A=(${K0C:-${LA[@]: -8}})\n"
       "  ( CA=(${K0A[@]}); NC=${#CA[@]}; plain25 0 $T0 k0w; [ -n \"$RK\" ] && echo \"0=$RK\" > $OUT/k0w.part ) &\n"
       "  P0=$!; log \"k0 WIRED ($(ngof 0) generic outputs): loader2 on ${K0A[*]} (exact path, seat C5)\"\n"
       "fi\n"
       "( CA=(${RA[@]}); NC=${#CA[@]}; right_chain ) &\nPR=$!\n")
out = src.replace(old, new)
md5 = hashlib.md5(src.encode()).hexdigest()[:8]
hdr = ("#!/bin/bash\n"
       "# chem_c5_port.sh: seat C5 (Nov 19 20:12Z), owner of the port chain since the lead's 20:00Z. A COPY of builder C4's\n"
       f"# tools/sh/chem_c4_port5.sh (md5 {md5}) with ONE change: the EXACT path (FAST=0, = port4's) builds a WIRED k0 in the\n"
       "# port (loader2 on the K0C cores) as the FAST path does, so k0-wired layouts (K0P57x, K0G10P57x, K0G11P57x) get an exact\n"
       "# port too. Written by tools/chem_c5_mkport.py. Everything else is port5's.\n")
assert out.startswith("#!/bin/bash\n")
out = hdr + out[len("#!/bin/bash\n"):]
open(os.path.join(here, "sh", "chem_c5_port.sh"), "w", newline="\n").write(out)
print(f"written sh/chem_c5_port.sh from chem_c4_port5.sh md5 {md5}")
