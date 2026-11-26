#!/usr/bin/env python3
# Audit k_47: relations/span.py's own steps (a)-(c) on its own rows WITHOUT the data line 5417 (imports span.py on purpose; not part of
# the independent verification).  The rows are those recorded from span.relation_rows by span_rerun.py (out/span_stage1_47.pkl).
# Without line 5417 the primes 6453089 and 184398551 occur in one line only (10574), so step (a) must remove them and that line.
# Expected if the audit's analysis is right: T has 9154 primes and the gcd of step (c) is 1 (no use of Iwasawa's theorem).
#   python3 -O span_reduced.py [LINE ...]        (default: 5417)
import sys, os, pickle, types, time
here=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(here,"..","..","relations"))
import span
from flint import fmpz_poly
t0=time.time(); p=47; drop=[int(a) for a in sys.argv[1:]] or [5417]
st=pickle.load(open(os.path.join(here,"out","span_stage1_47.pkl"),"rb"))
rows=[{q:fmpz_poly(v) for q,v in r.items()} for i,r in enumerate(st["rows"]) if i not in drop]
print("rows of span.py: %d; left out: %s; used: %d"%(len(st["rows"]),drop,len(rows)),flush=True)
T,g=span.prove(types.SimpleNamespace(N=p),rows,120,(4,400))          # the arguments of span.main
print("result: |T| = %d; 6453089 in T: %s; 184398551 in T: %s; gcd of step (c) = %s; %.0f s"%(len(T),6453089 in T,184398551 in T,g if g<10**9 else "%d digits"%len(str(g)),time.time()-t0))
open(os.path.join(here,"out","span_reduced_T.txt"),"w").write("".join("%d\n"%q for q in T))
