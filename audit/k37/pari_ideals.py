#!/usr/bin/env python3
# Audit k_37, another angle: the prime ideal factorizations and the Galois action are recomputed with PARI's number field functions
# (nfinit with the integral basis 1, eta_0..eta_35 written in powers of eta_0; idealprimedec; nfeltval; nfgaloisapply).
# No homomorphism to F_q is used here.  For a prime q, P_0 = the first prime ideal of idealprimedec, P_j = tau^(-j)(P_0) = s^j(P_0).
# For every row of the certificate that contains q:  V = (v_{P_j}(alpha))_j.  Checks: V has entries 0/1 only, sum = exponent of q,
# and {j : V_j = 1} = J_own + d  with ONE shift d for all rows (J_own from rows37.py).  That the shift does not depend on the row
# is the statement that the entries sum_{j in J} x^j of all rows refer to the same base prime ideal and the same generator s.
#   python3 pari_ideals.py q1 q2 ...      (default: 691 787 1013 1069 and some larger ones)
import sys, json, time
import cypari2
from field37 import Field, need, read_cert
pari=cypari2.Pari(); pari.allocatemem(4*10**9,silent=True)
p=37; F=Field(p); t0=time.time()
gal=json.load(open("out/galois.json")); den=pari(gal["den"]); x=pari("x")
f=pari.Pol([pari(c) for c in reversed(gal["f"])])
g=[pari.Pol([pari(c) for c in reversed(pl)])/den for pl in gal["polys"]]      # eta_i = g[i](eta_0)
B=[pari(1)]+g[:p-1]
nf=pari.nfinit([f,B]); need(nf[2]==pari(37)**72,"nfinit: discriminant"); print("nfinit: disc = 37^72; %.0f s"%(time.time()-t0),flush=True)
cert=read_cert("../../relations/k37_1/span_cert.txt"); rows=json.load(open("out/rows.json"))
allq=sorted(set(q for r in rows for q,J in r)); quiet=False
if len(sys.argv)>1 and sys.argv[1]=="range":          # python3 pari_ideals.py range A B : the primes number A..B-1 of the sorted list of all primes
    qs=allq[int(sys.argv[2]):int(sys.argv[3])]; quiet=True
else: qs=[int(a) for a in sys.argv[1:]] or [691,787,1013,1069,1669,36469,8418659,984681211]
byq={}
for n,row in enumerate(rows):
    for q,J in row: byq.setdefault(q,[]).append(n)
tot=0; budget=float(sys.argv[4]) if quiet and len(sys.argv)>4 else 1e9; done=0
for q in qs:
    if time.time()-t0>budget: break
    dec=pari.idealprimedec(nf,q); need(len(dec)==p and all(int(P[2])==1 and int(P[3])==1 for P in dec),"%d: not 37 primes of degree 1"%q)
    hn=[pari.idealhnf(nf,P) for P in dec]; P0=hn[0]; Pj=[]
    for j in range(p):
        I=pari.idealhnf(nf,pari.nfgaloisapply(nf,g[(-j)%p],P0)); k=[t for t in range(p) if hn[t]==I]; need(len(k)==1,"%d: conjugate not found"%q); Pj.append(dec[k[0]])
    need(len(set(str(pari.idealhnf(nf,P)) for P in Pj))==p,"%d: the conjugates are not distinct"%q)
    shifts=set(); nrow=0; kmax=0
    for n in byq[q]:
        (c,coef,fac),row=cert[n],rows[n]
        d=dict((a,b) for a,b in row)
        al=c+sum(a*g[i] for i,a in coef); V=[int(pari.nfeltval(nf,al,P)) for P in Pj]
        need(set(V)<={0,1} and sum(V)==dict(fac)[q],"%d: valuations %s"%(q,V))
        JV=set(j for j in range(p) if V[j]); ds=[s for s in range(p) if set((j+s)%p for j in d[q])==JV]
        need(ds,"%d: the pattern is not a shift of the audit's pattern"%q)
        shifts=shifts&set(ds) if nrow else set(ds); nrow+=1; kmax=max(kmax,sum(V))
        need(shifts,"%d: the shift depends on the row"%q)
    tot+=nrow; done+=1
    if not quiet: print("q = %d: %d rows; valuations 0/1 with the right number; patterns = audit's patterns shifted by d = %s for all rows (largest exponent %d); %.0f s"%(q,nrow,sorted(shifts),kmax,time.time()-t0),flush=True)
if quiet: print("primes number %s..%d of %d (q = %d..%d): %d primes, %d entries; all valuations 0/1 with the right number, and for every prime the patterns are the audit's patterns shifted by one d; %.0f s"%(sys.argv[2],int(sys.argv[2])+done-1,len(allq),qs[0],qs[done-1],done,tot,time.time()-t0),flush=True)
