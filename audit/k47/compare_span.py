#!/usr/bin/env python3
# Audit k_47: copied unchanged from audit/k41_k43 (general p).
# Audit k_41, k_43: the data recorded from relations/span.py (span_rerun.py) against the audit's own data.  Does not import span.py.
# After audit/k37/compare_span.py.
#  (1) periods: span's R_q satisfies the multiplication table modulo q (so it is a ring homomorphism O_K -> F_q), and it is a cyclic
#      shift R^span_i = R^own_(i+d_q) of the audit's homomorphism.
#  (2) rows: span's entry at q equals x^(-d_q) times the audit's entry (the same relation, another base prime ideal), for every row.
#  (3) step (c): for the matrices C recorded from span.py, Nm(det C) by nmdet2.nm_crt (other auxiliary primes, PARI, other bound) = D of span.py.
#   python3 compare_span.py p SECONDS       (part (3) is resumed from out/compare_blocks_P.json)
import sys, json, pickle, math, time, os
from fieldp import Field, need
from nmdet2 import nm_crt
sys.set_int_max_str_digits(1000000)
p=int(sys.argv[1]); budget=float(sys.argv[2]); F=Field(p); t0=time.time()
st=pickle.load(open("out/span_stage1_%d.pkl"%p,"rb")); own=json.load(open("out/periods_%d.json"%p)); rows_own=json.load(open("out/rows_%d.json"%p))
tab=[(i,j,F.prod[i][j][0],F.prod[i][j][1]) for i in range(p) for j in range(i,p)]
shift={}
for q,R in st["R"].items():
    need(len(R)==p and sum(R)%q==0 and all((R[i]*R[j]-const-sum(mk*Rk for mk,Rk in zip(m,R)))%q==0 for i,j,const,m in tab),"span's periods for %d are not a ring homomorphism"%q)
    Ro=own[str(q)]; ds=[d for d in range(p) if all(R[i]==Ro[(i+d)%p] for i in range(p))]
    need(len(ds)==1,"%d: span's periods are not a cyclic shift of the audit's"%q); shift[q]=ds[0]
need(set(shift)==set(int(q) for q in own),"the two sets of primes differ")
print("(1) periods of span.py: %d primes, all verified against the multiplication table, all cyclic shifts of the audit's homomorphisms; %.0f s"%(len(shift),time.time()-t0),flush=True)
red=lambda e:[e[k]-e[p-1] for k in range(p-1)]
need(len(st["rows"])==len(rows_own),"number of rows"); cnt=0
for rs,ro in zip(st["rows"],rows_own):
    need(set(rs)==set(q for q,J in ro),"support")
    for q,J in ro:
        e=[0]*p
        for j in J: e[(j-shift[q])%p]=1               # own J: phi_own(tau^j alpha) = 0;  phi_span = phi_own tau^d  =>  J_span = J_own - d
        v=list(rs[q])+[0]*(p-1-len(rs[q])); need(v==red(e),"entry at %d differs"%q); cnt+=1
print("(2) rows of span.py: %d rows, %d entries, all equal to the audit's entries up to the shift of the base prime ideal; %.0f s"%(len(rows_own),cnt,time.time()-t0),flush=True)
if not os.path.exists("out/span_blocks_%d.pkl"%p): raise SystemExit("(3) no recorded sets yet")
blocks=pickle.load(open("out/span_blocks_%d.pkl"%p,"rb")); fn="out/compare_blocks_%d.json"%p
res=json.load(open(fn)) if os.path.exists(fn) else {}
for b,B in enumerate(blocks):
    if str(b) in res or B["D"] is None: continue
    if time.time()-t0>budget: break
    n=B["n"]; fl=B["flat"]; C=[[fl[(i*n+j)*(p-1):(i*n+j+1)*(p-1)] for j in range(n)] for i in range(n)]
    D,k=nm_crt(C,p,flint_checks=1); need(abs(D)==B["D"],"set %d: D differs"%(b+1)); res[str(b)]=str(abs(D)); json.dump(res,open(fn,"w"))
    print("(3) set %d of span.py (n = %d): own Nm(det) %s (%d auxiliary primes); equal to the number of span.py; %.0f s"%(b+1,n,"= 0" if D==0 else "has %d digits"%len(str(abs(D))),k,time.time()-t0),flush=True)
g=0
for b in range(len(blocks)):
    if str(b) in res: g=math.gcd(g,int(res[str(b)])); print("    after set %d: gcd %s"%(b+1,"= %d"%g if g<10**9 else "has %d digits"%len(str(g))))
print("(3) sets compared: %d of %d; gcd of the compared ones %s"%(len(res),len(blocks),"= %d"%g if g<10**9 else "has %d digits"%len(str(g))))
