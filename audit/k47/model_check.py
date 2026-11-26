#!/usr/bin/env python3
# Audit k_41, k_43: the field model of fieldp.py against PARI; O_K = Z[1, eta_0..eta_(p-2)]; the prime above p; the thresholds.
# Adapted from audit/k37/model_check.py (constants 37, 1369, 72 replaced by p, p^2, 2(p-1); read and rewritten).
#   python3 model_check.py p
import sys, time
import cypari2
from fractions import Fraction
from fieldp import Field, need, is_prime_td
pari=cypari2.Pari(); pari.allocatemem(2*10**9,silent=True)
p=int(sys.argv[1]); F=Field(p,bruteforce=True); M=F.M; t0=time.time()
print("p = %d, M = %d, least primitive root g = %d, |H| = %d, g^p in H, H -> (Z/p)^* bijective, zeta -> zeta^g sends eta_i to eta_(i+1): checked in fieldp.py"%(p,M,F.g,len(F.H)))
print("H = <g^p>: checked (the definition of relations/verify.py);  E_i = g^i H")
E0=F.mul_matrix(0,[(0,1)])
f=pari.matrix(p,p,[x for r in E0 for x in r]).charpoly()
fs=pari("polsubcyclo(%d,%d)"%(M,p))
print("charpoly of multiplication by eta_0 == polsubcyclo(%d,%d): %s"%(M,p,f==fs)); need(f==fs,"charpoly")
need(pari.polisirreducible(f),"charpoly is reducible"); print("irreducible: True")
# Galois: the table is invariant under i -> i+1 (so tau: eta_i -> eta_(i+1) is a ring automorphism of the lattice)
ok=all(F.prod[(i+1)%p][(j+1)%p]==(F.prod[i][j][0],[F.prod[i][j][1][(k-1)%p] for k in range(p)]) for i in range(p) for j in range(p))
print("the table is invariant under the shift i -> i+1:",ok); need(ok,"shift")
# discriminant of the lattice Z[1, eta_0..eta_(p-2)]: det(Tr(b_i b_j)), Tr = trace of the regular representation
basis=[[1]+[0]*(p-1)]+[[0]*(1+k)+[1]+[0]*(p-2-k) for k in range(p-1)]
mats=[[[1 if a==b else 0 for b in range(p)] for a in range(p)]]+[F.mat[k] for k in range(p-1)]       # mats[t][k] = b_k * b_t
tr=[sum(mats[t][k][k] for k in range(p)) for t in range(p)]                                         # Tr(b_t)
G=[[sum(mats[j][i][t]*tr[t] for t in range(p)) for j in range(p)] for i in range(p)]                # Tr(b_i b_j)
dG=pari.matrix(p,p,[x for r in G for x in r]).matdet()
print("Tr(1) = %d, Tr(eta_k) = %s"%(tr[0],sorted(set(tr[1:]))))
print("discriminant of the lattice Z[1, eta_0, ..., eta_%d] = %s"%(p-2,pari.factor(dG)))
need(dG==pari(p)**(2*(p-1)),"lattice discriminant is not p^(2(p-1))")
print("conductor-discriminant formula: the p-1 nontrivial characters of Gal(K/Q) have conductor p^2, so d(K) = p^(2(p-1)) = lattice discriminant  =>  the lattice is O_K (index 1);  log rd = (2(p-1)/p) log p")
# PARI: maximal order at p only (nfbasis with the list [p]) and, if it finishes quickly, nfdisc
nd=pari("alarm(240, nfdisc(%s))"%f)                  # PARI's own alarm: returns an error object if the time is over
if nd.type()=="t_INT":
    print("PARI nfdisc(f) = %s; equal to p^(2(p-1)): %s"%(str(pari.factor(nd)).replace("\n",""),nd==pari(p)**(2*(p-1))))
    D=pari.poldisc(f); ind2=D//nd; ind=pari.sqrtint(ind2); need(ind*ind*nd==D,"index")
    print("index [O_K : Z[eta_0]] has %d digits; factors found by trial division up to 10^7: %s"%(len(str(ind)),str(pari.factor(ind,10**7)).replace("\n","")))
else: print("PARI nfdisc did not finish in 240 s")
# pi = prod_{h in H} (1 - zeta^h): N(pi) = Phi_M(1) = p
cnt=[0]*M; cnt[0]=1
for h in F.H:
    new=[0]*M
    for e in range(M):
        if cnt[e]: new[e]+=cnt[e]; new[(e+h)%M]-=cnt[e]
    cnt=new
mp=set(cnt[p*k] for k in range(1,p)); need(len(mp)==1,"pi: multiples of p"); const=cnt[0]-mp.pop(); m=[]
for k in range(p):
    s=set(cnt[e] for e in F.E[k]); need(len(s)==1,"pi: coset"); m.append(s.pop())
dpi=pari.matrix(p,p,[x for r in F.mul_matrix(const,[(k,m[k]) for k in range(p) if m[k]]) for x in r]).matdet()
print("pi = prod_{h in H}(1 - zeta^h) lies in Z[1, eta_i] (H-invariant coefficients); N(pi) = %s"%dpi); need(abs(int(dpi))==p,"norm of pi")
# threshold: min over primes l != p of l^(order of l mod p)
best=None; l=1
while best is None or l<best[0]:
    l+=1
    if l==p or not is_prime_td(l): continue
    fo=1
    while pow(l,fo,p)!=1: fo+=1
    if best is None or l**fo<best[0]: best=(l**fo,l,fo)
print("min over primes l != p of l^(order of l mod p) = %d (l = %d, order %d); searched all l < %d"%(best[0],best[1],best[2],best[0]))
small=[(l,next(fo for fo in range(1,p) if pow(l,fo,p)==1)) for l in range(2,best[0]) if is_prime_td(l) and l!=p]
print("orders modulo p of the primes below the threshold (l:order, l^order):"," ".join("%d:%d"%(l,fo) for l,fo in small))
print("  second smallest value of l^f:",sorted(l**fo for l,fo in [(l,next(fo for fo in range(1,p) if pow(l,fo,p)==1)) for l in range(2,3000) if is_prime_td(l) and l!=p])[:3])
print("%.0f s"%(time.time()-t0))
