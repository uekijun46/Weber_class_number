#!/usr/bin/env python3
# Audit k_41, k_43, step 1: every line of the certificate is checked again.  Nothing from relations/ is imported.
#   (i)   |N(alpha)| = |det of the regular representation of alpha| in the basis 1, eta_0..eta_(p-2), computed twice: FLINT fmpz_mat.det
#         and PARI matdet (both exact).  It must equal the product of the claimed factorization.  Every factor is prime (own deterministic
#         Miller-Rabin for numbers below 3*10^24, and PARI isprime) and is p or a split prime (q != p, q^(p-1) = 1 mod p^2).
#   (ii)  for every split prime q a ring homomorphism phi_q: O_K -> F_q.  It is FOUND as a common eigenvector: phi is a ring homomorphism
#         iff the column vector v = (phi(b_k))_k satisfies Mat(eta_i) v = phi(eta_i) v for all i (Mat(eta_i)[k] = coordinates of b_k eta_i).
#         Start with the eigenspace of Mat(eta_0) for a root rho_0 of the characteristic polynomial f of eta_0 modulo q; while its
#         dimension is > 1, intersect with an eigenspace of Mat(eta_1), Mat(eta_2), ...  (needed when q divides the index of Z[eta_0],
#         e.g. q = 19 < 43 for k_43, where O_K/q has no generator at all).  This is another construction than relations/span.py
#         (extension field F_q[t]/(m), zeta -> element of order M) and than audit/k37/rows37.py (powers of eta_0 + t eta_1).
#         However it was found, it is VERIFIED: R_i = phi(eta_i) satisfy sum_i R_i = 0 and R_i R_j = const_ij + sum_k m_ijk R_k (mod q)
#         for all i <= j (table of fieldp.py).  Since O_K = Z[1, eta_0..eta_(p-2)] (model_check.py), this proves that eta_i -> R_i
#         defines a ring homomorphism O_K -> F_q.  The p cyclic shifts of R are distinct (checked): they are the p homomorphisms.
#   (iii) J = {j : c + sum_i a_i R_(i+j) = 0 mod q};  required: |J| = exponent of q in the norm.
# Output: out/rows_P.json (list of rows; a row is a list of [q, J]),  out/periods_P.json (q -> R),  text on stdout.
#   python3 rows.py p
import sys, json, time, collections
import cypari2
from flint import fmpz_mat, nmod_poly, nmod_mat
from fieldp import Field, need, read_cert, is_prime_mr, CERT
pari=cypari2.Pari(); pari.allocatemem(10**9,silent=True)
p=int(sys.argv[1]); F=Field(p); t0=time.time()
cert=read_cert(CERT(p),p); print("field k_%d; lines: %d"%(p,len(cert)),flush=True)
# (i)
primes=set(); maxabs=0; nz=collections.Counter(); cs=collections.Counter(); coefvals=collections.Counter(); big=collections.Counter(); occ=collections.Counter(); seen=set(); withp=0
for n,(c,coef,fac) in enumerate(cert):
    need(all(0<=i<p and a!=0 for i,a in coef) and len(set(i for i,a in coef))==len(coef),"coefficients, line %d"%n)
    need(len(set(r for r,e in fac))==len(fac) and all(e>=1 for r,e in fac),"factorization, line %d"%n)
    prod=1
    for r,e in fac:
        need(is_prime_mr(r) and bool(pari.isprime(r)),"%d is not prime, line %d"%(r,n)); need(r==p or F.split(r),"%d is not p or a split prime, line %d"%(r,n)); prod*=r**e
        if r!=p: primes.add(r); occ[r]+=1
    A=F.mul_matrix(c,coef); flat=[x for row in A for x in row]
    d1=int(fmpz_mat(p,p,flat).det()); d2=int(pari.matrix(p,p,flat).matdet())
    need(d1==d2,"the two determinants differ, line %d"%n); need(abs(d1)==prod,"norm, line %d: %d vs %d"%(n,d1,prod))
    maxabs=max(maxabs,abs(d1)); nz[len(coef)]+=1; cs[c]+=1; withp+=any(r==p for r,e in fac)
    for i,a in coef: coefvals[a]+=1
    big[sum(1 for r,e in fac if r>=10**6)]+=1
    key=(c,tuple(sorted(coef))); need(key not in seen,"duplicate line %d"%n); seen.add(key)
print("(i) norms (two exact determinants) and factorizations: %d lines, all equal; split primes: %d; largest norm has %d digits; %.0f s"%(len(cert),len(primes),len(str(maxabs)),time.time()-t0),flush=True)
print("    shape: number of nonzero coefficients %s; c %s; coefficient values %s; number of prime factors >= 10^6 per line %s; lines whose norm is divisible by p: %d"%(sorted(nz.items()),sorted(cs.items()),sorted(coefvals.items()),sorted(big.items()),withp))
print("    largest prime %d; primes below 10^6: %d; smallest number of lines in which a prime occurs: %d; no duplicate lines"%(max(primes),sum(1 for q in primes if q<10**6),min(occ.values())))
# (ii)
E0=F.mat[0]; fpol=[int(x) for x in pari.Vecrev(pari.matrix(p,p,[x for r in E0 for x in r]).charpoly())]
tab=[(i,j,F.prod[i][j][0],F.prod[i][j][1]) for i in range(p) for j in range(i,p)]
flatmat=[[x for r in F.mat[i] for x in r] for i in range(p)]
def hom(q):
    roots=sorted(int(r) for r,e in nmod_poly(fpol,q).roots())
    need(roots,"%d: f has no root"%q)
    W=None; steps=0
    for i in range(p):
        Mi=nmod_mat(p,p,flatmat[i],q); found=False
        for rho in roots:
            D=Mi-nmod_mat(p,p,[rho if a==b else 0 for a in range(p) for b in range(p)],q)
            if W is not None: D=D*W
            X,nul=D.nullspace()
            if nul>0:
                X=nmod_mat(X.nrows(),nul,[int(X[a,b]) for a in range(X.nrows()) for b in range(nul)],q)
                W=X if W is None else W*X; found=True; break
        need(found,"%d: no common eigenvector"%q); steps+=1
        if W.ncols()==1: break
    need(W.ncols()==1,"%d: the common eigenspace has dimension %d"%(q,W.ncols()))
    v=[int(W[k,0]) for k in range(p)]; need(v[0]%q!=0,"%d: phi(1) = 0"%q); iv=pow(v[0],-1,q)
    R=[v[1+k]*iv%q for k in range(p-1)]; R.append(-sum(R)%q)
    return R,steps
def check_hom(q,R):
    if len(R)!=p or sum(R)%q: return False
    for i,j,const,m in tab:
        if (R[i]*R[j]-const-sum(mk*Rk for mk,Rk in zip(m,R)))%q: return False
    return len(set(tuple(R[j:]+R[:j]) for j in range(p)))==p
R={}; st=collections.Counter(); nvals=collections.Counter()
for n,q in enumerate(sorted(primes)):
    R[q],steps=hom(q); st[steps]+=1
    need(check_hom(q,R[q]),"%d: not a ring homomorphism"%q)
    if len(set(R[q]))<p: nvals[q]=len(set(R[q]))
    if n%1000==0: print("  prime %d of %d, %.0f s"%(n,len(primes),time.time()-t0),flush=True)
print("(ii) ring homomorphisms O_K -> F_q verified by the multiplication table for %d primes (number of eigenspace steps: %s); %.0f s"%(len(R),sorted(st.items()),time.time()-t0),flush=True)
print("    primes where the values R_0..R_(p-1) are not distinct (q divides the index of Z[eta_0]): %d, the smallest: %s"%(len(nvals),sorted(nvals.items())[:12]))
# (iii)
rows=[]; ksum=collections.Counter(); nent=0
for n,(c,coef,fac) in enumerate(cert):
    row=[]
    for q,k in fac:
        if q==p: continue
        Rq=R[q]; J=[j for j in range(p) if (c+sum(a*Rq[(i+j)%p] for i,a in coef))%q==0]
        need(len(J)==k,"line %d, prime %d: %d prime ideals contain the element, exponent %d"%(n,q,len(J),k))
        ksum[k]+=1; nent+=1; row.append([q,J])
    rows.append(row)
print("(iii) every exponent equals the number of prime ideals containing the element: %d lines, %d entries; exponents: %s; entries with all p prime ideals: %d; %.0f s"%(len(rows),nent,sorted(ksum.items()),ksum.get(p,0),time.time()-t0))
json.dump(rows,open("out/rows_%d.json"%p,"w")); json.dump({str(q):v for q,v in R.items()},open("out/periods_%d.json"%p,"w"))
