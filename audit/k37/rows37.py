#!/usr/bin/env python3
# Audit k_37, step 1: every line of the certificate is checked again, without relations/span.py and relations/verify.py.
#   (i)  |N(alpha)| = |det of the regular representation of alpha| in the basis 1, eta_0..eta_35 (PARI matdet, exact) = product of the
#        claimed factorization; every factor is prime (PARI isprime, a proof) and is 37 or a split prime (q^36 = 1 mod 1369).
#   (ii) for every split prime q: a ring homomorphism phi_q: O_K -> F_q.  O_K = Z[1, eta_0..eta_35] (index 1, model_check.py).
#        phi_q is found by linear algebra over F_q (theta = eta_0 + t eta_1 generates O_K/q, a root rho of its characteristic polynomial,
#        eta_i as polynomials in theta modulo q).  However it was found, it is VERIFIED: R_i = phi_q(eta_i) satisfy sum R_i = 0 and
#        R_i R_j = const_ij + sum_k m_ijk R_k (mod q) for all i <= j.  This proves that eta_i -> R_i is a ring homomorphism.
#   (iii) J = {j : c + sum a_i R_(i+j) = 0 mod q};  required: |J| = exponent of q in the norm.
# Output: out/rows.json  (list of rows; a row is a list of [q, J]),  out/periods.json (q -> R).
import sys, json, time
import cypari2
from field37 import Field, need, read_cert
pari=cypari2.Pari(); pari.allocatemem(10**9)
p=37; F=Field(p); t0=time.time()
cert=read_cert(sys.argv[1]); print("lines:",len(cert),flush=True)
basis=[[1]+[0]*(p-1)]+[[0]*(1+k)+[1]+[0]*(p-2-k) for k in range(p-1)]
def mul(u,v):
    Mv=F.mul_matrix(v); return [sum(u[k]*Mv[k][t] for k in range(p)) for t in range(p)]
# (i)
primes=set(); maxabs=0
for n,(c,coef,fac) in enumerate(cert):
    need(all(0<=i<p and a!=0 for i,a in coef) and len(set(i for i,a in coef))==len(coef),"coefficients, line %d"%n)
    need(len(set(r for r,e in fac))==len(fac) and all(e>=1 for r,e in fac),"factorization, line %d"%n)
    prod=1
    for r,e in fac:
        need(bool(pari.isprime(r)),"%d is not prime, line %d"%(r,n)); need(r==p or F.split(r),"%d is not 37 or split, line %d"%(r,n)); prod*=r**e
        if r!=p: primes.add(r)
    d=int(pari.matrix(p,p,[x for row in F.mul_matrix(F.elem(c,coef)) for x in row]).matdet())
    need(abs(d)==prod,"norm, line %d: %d vs %d"%(n,d,prod))
    maxabs=max(maxabs,abs(d))
print("(i) norms and factorizations: %d lines, all equal; primes other than 37: %d; largest norm has %d digits; %.0f s"%(len(cert),len(primes),len(str(maxabs)),time.time()-t0),flush=True)
# (ii)
tab=[(i,j,F.prod[i][j]) for i in range(p) for j in range(i,p)]
def hom(q):
    for t in range(0,50):
        th=F.elem(0,[(0,1)]+([(1,t)] if t else []))
        pw=[basis[0]]
        for k in range(1,p): pw.append([x%q for x in mul(pw[-1],th)])
        P=pari.matrix(p,p,[x for r in pw for x in r])*pari.Mod(1,q)
        if pari.matdet(P)==0: continue                # q divides the index of Z[theta]
        Pinv=P**(-1)
        pw.append([x%q for x in mul(pw[-1],th)])      # theta^p = sum_k w_k theta^k: w = (theta^p) * Pinv
        w=pari.vector(p,pw[p])*Pinv
        f=pari('x')**p-pari.Pol(list(reversed([w[k] for k in range(p)])))
        roots=pari.polrootsmod(pari.lift(f),q)
        need(len(roots)==p,"%d: the polynomial of theta does not have 37 roots"%q)
        rho=int(pari.lift(roots[0])); rp=[pow(rho,k,q) for k in range(p)]
        R=[sum(int(pari.lift(Pinv[1+i,k]))*rp[k] for k in range(p))%q for i in range(p-1)]
        R.append(-sum(R)%q)
        return R,t
    need(False,"%d: no generator found"%q)
def check_hom(q,R):
    if sum(R)%q: return False
    for i,j,(const,m) in tab:
        if (R[i]*R[j]-const-sum(m[k]*R[k] for k in range(p)))%q: return False
    return len(set(tuple(R[j:]+R[:j]) for j in range(p)))==p
R={}; tused={}
for n,q in enumerate(sorted(primes)):
    R[q],t=hom(q); tused[t]=tused.get(t,0)+1
    need(check_hom(q,R[q]),"%d: not a ring homomorphism"%q)
    if n%500==0: print("  prime %d of %d, %.0f s"%(n,len(primes),time.time()-t0),flush=True)
print("(ii) ring homomorphisms O_K -> F_q verified for %d primes (generator eta_0 + t eta_1, t: %s); %.0f s"%(len(R),sorted(tused.items()),time.time()-t0),flush=True)
# (iii)
rows=[]; full=0; ksum={}
for n,(c,coef,fac) in enumerate(cert):
    row=[]
    for q,k in fac:
        if q==p: continue
        Rq=R[q]; J=[j for j in range(p) if (c+sum(a*Rq[(i+j)%p] for i,a in coef))%q==0]
        need(len(J)==k,"line %d, prime %d: %d prime ideals contain the element, exponent %d"%(n,q,len(J),k))
        ksum[k]=ksum.get(k,0)+1
        if k==p: full+=1
        row.append([q,J])
    rows.append(row)
print("(iii) every exponent equals the number of prime ideals containing the element: %d lines; exponents: %s; entries with all 37 prime ideals: %d; %.0f s"%(len(rows),sorted(ksum.items()),full,time.time()-t0))
json.dump(rows,open("out/rows.json","w")); json.dump({str(q):v for q,v in R.items()},open("out/periods.json","w"))
