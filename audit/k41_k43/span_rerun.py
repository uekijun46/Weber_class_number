#!/usr/bin/env python3
# Audit k_41, k_43: REPRODUCTION of relations/span.py, with hooks that record intermediate data (after audit/k37/span_staged.py).
# This script imports span.py on purpose; it is NOT part of the independent verification (rows.py, elim.py, core.py).
# span.main() is called unchanged with the arguments  CERT --primes out/span_T_rerun_P.txt .  Hooks:
#   span.periods is wrapped (records R_q);  span.relation_rows is wrapped (records the rows);
#   span.nmod_mat is wrapped (records the coefficient array of each set of step (c));  span.math is replaced by an object whose gcd records D.
# The computations of span.py are not changed.  Output: out/span_stage1_P.pkl (rows and periods), out/span_blocks_P.pkl, text on stdout.
#   python3 -O span_rerun.py p
import sys, os, pickle, time, math, types
here=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(here,"..","..","relations"))
import span
sys.set_int_max_str_digits(1000000)
p=int(sys.argv[1]); t0=time.time()
cert=os.path.join(here,"..","..","relations","k%d_1"%p,"span_cert.txt")
rec={}; orig_periods=span.periods
def periods(F_,q): rec[q]=orig_periods(F_,q); return rec[q]
span.periods=periods
orig_rr=span.relation_rows
def relation_rows(F_,rels):
    rows=orig_rr(F_,rels)
    pickle.dump({"rows":[{q:[int(c) for c in v.coeffs()] for q,v in r.items()} for r in rows],"R":rec},open(os.path.join(here,"out","span_stage1_%d.pkl"%p),"wb"))
    print("[audit hook] rows and periods recorded; %.0f s"%(time.time()-t0),flush=True); return rows
span.relation_rows=relation_rows
blocks=[]; orig_mat=span.nmod_mat; last=[None]
def nmod_mat(*a):
    if len(a)==4 and a[1]==p-1 and a[0]>p-1 and a[2] is not last[0]:
        last[0]=a[2]; blocks.append({"flat":list(a[2]),"n":math.isqrt(a[0]),"D":None})
    return orig_mat(*a)
span.nmod_mat=nmod_mat
def gcd(g,D):
    blocks[-1]["D"]=D; pickle.dump(blocks,open(os.path.join(here,"out","span_blocks_%d.pkl"%p),"wb")); return math.gcd(g,D)
span.math=types.SimpleNamespace(prod=math.prod,isqrt=math.isqrt,gcd=gcd)
sys.argv=["span.py",cert,"--primes",os.path.join(here,"out","span_T_rerun_%d.txt"%p)]
span.main()
print("[audit hook] total %.0f s; sets recorded: %d"%(time.time()-t0,len(blocks)))
