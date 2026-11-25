#!/usr/bin/env python3
# Audit k_37, corroboration only (not a proof): the rank of the core left by elim_own.py modulo some maximal ideals of A = Z[x]/(Phi_37).
# L' = A^n holds iff the rank is n modulo EVERY maximal ideal; here only a few are tested:
#   (37, x - 1): residue field F_37;   (l, x - w^k) for l = 149, 223, 593 and all k = 1..36: residue field F_l;   (2): residue field F_2^36.
import sys, pickle, time
import cypari2
pari=cypari2.Pari(); pari.allocatemem(2*10**9,silent=True)
p=37; st=pickle.load(open("out/elim_state.pkl","rb")); cols=st["cols"]; n=len(cols); alive=st["alive"]; t0=time.time()
pad=lambda v:list(v)+[0]*(p-1-len(v))
R=[[pad(st["rows"][i].get(q,[])) for q in cols] for i in alive]
M=pari.matrix(len(R),n,[sum(e)%37 for r in R for e in r])*pari.Mod(1,37)
print("(37, x-1): rank %d of %d (%d rows)"%(int(pari.matrank(M)),n,len(R)),flush=True)
for l in (149,223,593):
    w=next(pow(a,(l-1)//p,l) for a in range(2,l) if pow(a,(l-1)//p,l)!=1); rk=[]
    for k in range(1,p):
        pw=[pow(w,k*d,l) for d in range(p-1)]
        M=pari.matrix(len(R),n,[sum(c*pw[d] for d,c in enumerate(e))%l for r in R for e in r])*pari.Mod(1,l); rk.append(int(pari.matrank(M)))
    print("(%d, x - w^k), k = 1..36: ranks %s (%.0f s)"%(l,sorted(set(rk)),time.time()-t0),flush=True)
Phi2=pari.polcyclo(37)*pari.Mod(1,2); need=pari.polisirreducible(Phi2); print("Phi_37 irreducible modulo 2:",bool(need))
sub=R[:400]; x=pari('x')
M=pari.matrix(len(sub),n,[pari.Mod(pari.Pol(list(reversed(e)))*pari.Mod(1,2),Phi2) for r in sub for e in r])
print("(2): rank %d of %d over F_2^36 (first %d rows; %.0f s)"%(int(pari.matrank(M)),n,len(sub),time.time()-t0))
