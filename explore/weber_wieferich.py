# "Wieferich-type" necessary condition for l | h (cyclotomic Z_2-extension), l = 1 mod 2^r:
#   l | h_psi  ==>  l | B_{1, psi*omega^{-1}}  <==>  sum_{a < M*l, a odd, l does not divide a} psi(a) q_l(a) = 0 (mod l),
# where psi is an even character of order 2^r and conductor M = 2^(r+2), and q_l(a) = (a^(l-1)-1)/l is the Fermat quotient.
import numpy as np, math, sys
from sympy import primerange, primitive_root
def powmod(a,e,m):
    r=np.ones_like(a); b=a%m
    while e:
        if e&1: r=(r*b)%m
        b=(b*b)%m; e>>=1
    return r
def run(r,L):
    M=2**(r+2); n=2**r
    idx={}                                   # a = +-5^k mod M  ->  k
    x=1
    for k in range(n):
        idx[x]=k; idx[M-x]=k; x=(x*5)%M
    kk=np.array([idx.get(u,-1) for u in range(M)])
    zeros=[]; expect=0.0; cnt=0
    for l in primerange(3,L):
        if l%n!=1: continue
        cnt+=1; expect+=(n//2)/l
        a=np.arange(1,M*l,2,dtype=np.int64); a=a[a%l!=0]
        q=((powmod(a,l-1,l*l)-1)//l)%l
        T=np.zeros(n,dtype=np.int64); np.add.at(T,kk[a%M],q); T%=l
        g=primitive_root(l); z=pow(g,(l-1)//n,l)
        for j in range(1,n,2):
            w=sum(int(T[k])*pow(z,j*k,l) for k in range(n))%l
            if w==0: zeros.append((l,j))
    print(f"r={r}: primes l = 1 mod {n} below {L}: {cnt}; expected number of vanishings ~ {expect:.2f}; found {len(zeros)}: {zeros}",flush=True)
for r,L in ((2,50000),(3,50000),(4,50000),(5,50000)):
    run(r,L)
