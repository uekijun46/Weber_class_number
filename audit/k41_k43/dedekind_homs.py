#!/usr/bin/env python3
# Audit k_41, k_43, third construction of the homomorphisms (for the primes q that do not divide the index [O_K : Z[eta_0]]):
#   eta_i = g_i(eta_0) with g_i in (1/d) Z[x]  (computed exactly from the multiplication table; checked numerically, mpmath 3000 digits,
#   against the definition eta_i = sum_{e in E_i} exp(2 pi i e / p^2)).  For q not dividing d and a root rho of f = charpoly(eta_0) modulo q,
#   eta_0 -> rho defines the ring homomorphism Z[eta_0] -> F_q, which extends uniquely to O_K (contained in Z[1/d][eta_0]): eta_i -> g_i(rho).
#   Required: f has p distinct roots modulo q, and the p vectors (g_i(rho))_i are exactly the p cyclic shifts of the audit's R_q (rows.py).
# Nothing from relations/ is imported.     python3 dedekind_homs.py p
import sys, json, time
import cypari2, mpmath as mp
from flint import nmod_poly
from fieldp import Field, need
pari=cypari2.Pari(); pari.allocatemem(2*10**9,silent=True)
p=int(sys.argv[1]); F=Field(p); M=F.M; t0=time.time()
pw=[[1]+[0]*(p-1)]; E0=F.mat[0]
for k in range(1,p): pw.append([sum(pw[-1][s]*E0[s][t] for s in range(p)) for t in range(p)])
P=pari.matrix(p,p,[x for r in pw for x in r]); need(P.matdet()!=0,"powers of eta_0"); Pinv=P**(-1); den=pari.denominator(Pinv); Pi=Pinv*den; deni=int(den)
polys=[]
for i in range(p):
    v=F.elem(0,[(i,1)]); polys.append([int(sum(v[t]*Pi[t,k] for t in range(p))) for k in range(p)])
mp.mp.dps=3000
num=[sum(mp.cos(2*mp.pi*e/M) for e in F.E[i]) for i in range(p)]; x0=num[0]; xp=[mp.mpf(1)]
for k in range(1,p): xp.append(xp[-1]*x0)
err=max(abs(sum(polys[i][k]*xp[k] for k in range(p))/deni-num[i]) for i in range(p)); need(err<mp.mpf(10)**(-1000),"g_i(eta_0) != eta_i numerically")
fpol=[int(c) for c in pari.Vecrev(pari.matrix(p,p,[c for r in E0 for c in r]).charpoly())]
ind=abs(int(P.matdet())); D=int(pari.poldisc(pari.Pol(list(reversed(fpol))))); need(D==ind*ind*p**(2*(p-1)) and ind%deni==0,"disc(f) != index^2 p^(2(p-1))")
print("k_%d: index of Z[eta_0] in O_K = |det(powers of eta_0)| has %d digits, disc(f) = index^2 * %d^%d; the common denominator d of the g_i has %d digits and divides the index (same prime divisors: exponent and order of the finite group O_K/Z[eta_0]); numerical check of g_i: largest error %s; %.0f s"%(p,len(str(ind)),p,2*(p-1),len(str(deni)),mp.nstr(err,3),time.time()-t0),flush=True)
Rown=json.load(open("out/periods_%d.json"%p)); allq=sorted(int(q) for q in Rown); skipped=[]; done=0
for q in allq:
    if ind%q==0: skipped.append(q); continue
    R=Rown[str(q)]; shiftof={tuple(R[j:]+R[:j]):j for j in range(p)}
    roots=[int(r) for r,e in nmod_poly(fpol,q).roots()]; need(len(set(roots))==p,"%d: f does not have p distinct roots"%q)
    pm=[[c%q for c in pl] for pl in polys]; di=pow(deni,-1,q); js=[]
    for rho in roots:
        vec=[]
        for pl in pm:
            a=0
            for c in reversed(pl): a=(a*rho+c)%q
            vec.append(a*di%q)
        need(vec[0]==rho and tuple(vec) in shiftof,"%d: g_i(rho) is not a cyclic shift of the audit's homomorphism"%q); js.append(shiftof[tuple(vec)])
    need(sorted(js)==list(range(p)),"%d: the roots do not give the p shifts"%q); done+=1
print("k_%d: %d primes not dividing the index: the homomorphisms eta_i -> g_i(rho) are exactly the p cyclic shifts of the audit's homomorphism; %.0f s"%(p,done,time.time()-t0))
print("primes of T dividing the index (left to pari_homs.py): %d: %s"%(len(skipped),skipped))
json.dump(skipped,open("out/index_primes_%d.json"%p,"w"))
