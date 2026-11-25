#!/usr/bin/env python3
# Audit k_37: comparison of the data recorded from span.py (span_staged.py) with the audit's own data.  Does not import span.py.
#  (1) periods: span's R_q satisfies the multiplication table modulo q (so it is a ring homomorphism O_K -> F_q), and it is a cyclic
#      shift R^span_i = R^own_(i+d_q) of the audit's homomorphism.
#  (2) rows: span's entry at q equals x^(-d_q) times the audit's entry (the same relation, another base prime ideal), for every row.
#  (3) step (c): for the matrices C recorded from span.py, Nm(det C) by nmdet.nm_crt (other auxiliary primes, PARI, other bound) = D of span.py.
import sys, json, pickle, math, time
from field37 import Field, need
from nmdet import nm_crt, zmatrix, hadamard_bound
sys.set_int_max_str_digits(1000000)
p=37; F=Field(p); t0=time.time()
st=pickle.load(open("out/span_stage1.pkl","rb")); own=json.load(open("out/periods.json")); rows_own=json.load(open("out/rows.json"))
tab=[(i,j,F.prod[i][j]) for i in range(p) for j in range(i,p)]
shift={}
for q,R in st["R"].items():
    need(sum(R)%q==0 and all((R[i]*R[j]-const-sum(m[k]*R[k] for k in range(p)))%q==0 for i,j,(const,m) in tab),"span's periods for %d are not a ring homomorphism"%q)
    Ro=own[str(q)]; ds=[d for d in range(p) if all(R[i]==Ro[(i+d)%p] for i in range(p))]
    need(len(ds)==1,"%d: span's periods are not a cyclic shift of the audit's"%q); shift[q]=ds[0]
print("(1) periods of span.py: %d primes, all verified against the multiplication table, all cyclic shifts of the audit's homomorphisms"%len(shift))
def red(e):                                           # list of length p (coefficients of x^0..x^(p-1)) -> length p-1 modulo Phi_p
    return [e[k]-e[p-1] for k in range(p-1)]
need(len(st["rows"])==len(rows_own),"number of rows"); cnt=0
for rs,ro in zip(st["rows"],rows_own):
    need(set(rs)==set(q for q,J in ro),"support")
    for q,J in ro:
        e=[0]*p
        for j in J: e[(j-shift[q])%p]=1               # own J: phi_own(tau^j alpha)=0;  phi_span = phi_own tau^d  =>  J_span = J_own - d
        v=list(rs[q])+[0]*(p-1-len(rs[q])); need(v==red(e),"entry at %d differs"%q); cnt+=1
print("(2) rows of span.py: %d rows, %d entries, all equal to the audit's entries up to the shift of the base prime ideal"%(len(rows_own),cnt))
blocks=pickle.load(open("out/span_blocks.pkl","rb")); g=0
for b,B in enumerate(blocks):
    n=B["n"]; fl=B["flat"]; C=[[fl[(i*n+j)*(p-1):(i*n+j+1)*(p-1)] for j in range(n)] for i in range(n)]
    D,k=nm_crt(C,p); g=math.gcd(g,abs(D))
    print("(3) set %d of span.py (n = %d): own Nm(det) %s; span.py: %s; equal up to sign: %s; gcd so far %s; %.0f s"%(b+1,n,"= 0" if D==0 else "has %d digits"%len(str(abs(D))),"= 0" if B["D"]==0 else "has %d digits"%len(str(B["D"])),abs(D)==B["D"],"= %d"%g if g<10**9 else "has %d digits"%len(str(g)),time.time()-t0),flush=True)
    need(abs(D)==B["D"],"D differs")
need(g==1,"gcd")
print("(3) gcd = 1 for the matrices of span.py")
