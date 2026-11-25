#!/usr/bin/env python3
# Audit k_37: relations/span.py (imported on purpose: the object under test) on modified inputs, starting from the rows recorded by
# span_staged.py stage 1 (out/span_stage1.pkl).  Modes:
#   core130      : prove(F, rows, 130, (4,400))   (another stopping point of step (b))
#   nocap        : prove(F, rows, 120, None)      (no row is dropped)
#   half         : only the first 3000 rows (fewer rows than primes)  -> must NOT pass
#   one691       : all rows that contain 691 are removed except the first one -> 691 occurs once, step (a) must remove it from T
#   no691        : all rows that contain 691 are removed -> 691 is not in T
#   zerorow K    : K rows are replaced by copies of other rows (duplicates)  -> harmless, must still pass or fail safely
import sys, os, pickle, time
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","relations"))
import span
from flint import fmpz_poly
mode=sys.argv[1]; t0=time.time(); F=span.K(37,1)
st=pickle.load(open("out/span_stage1.pkl","rb")); rows=[{q:fmpz_poly(v) for q,v in r.items()} for r in st["rows"]]
core,cap=120,(4,400)
if mode=="core130": core=130
elif mode=="nocap": cap=None
elif mode=="half": rows=rows[:3000]
elif mode=="one691":
    idx=[i for i,r in enumerate(rows) if 691 in r]; rows=[r for i,r in enumerate(rows) if i not in set(idx[1:])]; print("rows with 691: %d; kept 1"%len(idx))
elif mode=="no691":
    rows=[r for r in rows if 691 not in r]
elif mode=="dup":
    k=int(sys.argv[2]); rows=rows[:len(rows)-k]+[dict(r) for r in rows[:k]]
print("mode %s: %d rows"%(mode,len(rows)),flush=True)
try:
    T,g=span.prove(F,rows,core,cap)
    print("RESULT mode %s: gcd = %s; |T| = %d; 691 in T: %s; smallest %s; %.0f s"%(mode,g if g<10**9 else "(%d digits)"%len(str(g)),len(T),691 in T,T[:4],time.time()-t0))
except SystemExit as e:
    print("RESULT mode %s: %s; %.0f s"%(mode,e,time.time()-t0))
