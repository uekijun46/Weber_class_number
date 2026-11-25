#!/usr/bin/env python3
# Audit k_37: REPRODUCTION of relations/span.py in two stages (each below 10 minutes), with hooks that record intermediate data.
# This script imports span.py on purpose; it is NOT part of the independent verification (rows37.py, elim_own.py, core_own.py).
#   stage 1: span.relation_rows(span.read_cert(...))  -> out/span_stage1.pkl  (rows as coefficient lists, periods R_q of span.periods)
#   stage 2: span.prove(F, rows, 120, (4,400)) as in span.main  -> out/span_blocks.pkl  (the matrices C of step (c) and the numbers D)
#   stage 3: the lines "UNDER GRH" as in span.main
# Hooks: span.periods is wrapped (records R_q); span.nmod_mat is wrapped (records the coefficient array of each block);
# span.math is replaced by a module-like object whose gcd records D.  The computations of span.py are not changed.
import sys, os, pickle, time, math, types
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","relations"))
import span
from flint import fmpz_poly
sys.set_int_max_str_digits(1000000)
cert=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","relations","k37_1","span_cert.txt")
stage=sys.argv[1]; t0=time.time(); F=span.K(37,1)
if stage=="1":
    rec={}; orig=span.periods
    def periods(F_,q): rec[q]=orig(F_,q); return rec[q]
    span.periods=periods
    rows=span.relation_rows(F,span.read_cert(F,cert))
    pickle.dump({"rows":[{q:[int(c) for c in v.coeffs()] for q,v in r.items()} for r in rows],"R":rec},open("out/span_stage1.pkl","wb"))
    print("stage 1: %.0f s"%(time.time()-t0))
elif stage=="2":
    st=pickle.load(open("out/span_stage1.pkl","rb")); rows=[{q:fmpz_poly(v) for q,v in r.items()} for r in st["rows"]]
    blocks=[]; orig_mat=span.nmod_mat
    def nmod_mat(*a):
        if len(a)==4 and a[1]==36 and a[0]>36 and (not blocks or blocks[-1]["flat"]!=a[2]): blocks.append({"flat":list(a[2]),"n":math.isqrt(a[0]),"D":None})
        return orig_mat(*a)
    span.nmod_mat=nmod_mat
    def gcd(g,D):
        blocks[-1]["D"]=D; pickle.dump(blocks,open("out/span_blocks.pkl","wb")); return math.gcd(g,D)
    span.math=types.SimpleNamespace(prod=math.prod,isqrt=math.isqrt,gcd=gcd)
    T,g=span.prove(F,rows,120,(4,400))
    print("gcd = %d; |T| = %d; stage 2: %.0f s"%(g,len(T),time.time()-t0))
    open("out/span_T_rerun.txt","w").write("".join("%d\n"%q for q in T))
elif stage=="3":
    from bound_grh import bound_grh
    from flint import arb
    T=[int(l) for l in open("out/span_T_rerun.txt")]; hm=F.hmin(); best=None
    for c2 in range(8,29):
        c=arb(c2)/2; G_up,S_low,B,hb=bound_grh(F,T,c)
        print("  c=%s: G <= %s, sum >= %s, B >= %s, h < %s"%(c2/2,span.up(G_up,10),span.down(S_low,10),span.down(B,7),"-" if hb is None else span.up(hb,4)))
        if hb is not None and (best is None or hb<best[0]): best=(hb,c2/2,G_up,S_low,B)
    hb,c,G_up,S_low,B=best
    print("UNDER GRH (Miller, Lemma 5.2): c=%s: G(c) <= %s, prime sum >= %s, B >= %s, h < %s; least possible class number > 1: %d"%(c,span.up(G_up,10),span.down(S_low,10),span.down(B,7),span.up(hb,4),hm))
