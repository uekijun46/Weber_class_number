# Exploration: a faster version of the per-prime search of relations/search.py (untrusted search step).
# Changes: (1) before the full factorization, the split primes below a small bound and the powers of p are removed from the
# norm divided by q; if the rest is a prime number >= q, the candidate is rejected without factoring; (2) LLL with a smaller
# parameter delta (optional).
import sys, time, math, flint, cypari2
sys.path.insert(0,"../../relations"); sys.path.insert(0,".")
from field import Field
def setup(p,n,Y=200000):
    F=Field(p,n); F.fixed()
    small=[r for r in range(2,Y) if F.split(r) and flint.fmpz(r).is_prime()]
    P=1
    for r in small: P*=r
    return F,P
def targeted2(F,P,pari,q,delta=0.75,maxcand=400,maxl1=16):
    N=F.N; R=F.periods_mod(q); B=flint.fmpz_mat(N,N); B[0,0]=q
    for i in range(N-1): B[i+1,0]=(-R[i])%q; B[i+1,i+1]=1
    L=B.lll(delta=delta,eta=0.51); vs=[[int(L[i,j]) for j in range(N)] for i in range(N)]
    vs.sort(key=lambda v:sum(abs(x) for x in v)); short=vs[:14]; cands=list(vs)
    for i in range(len(short)):
        for j in range(i+1,len(short)):
            cands.append([a+b for a,b in zip(short[i],short[j])]); cands.append([a-b for a,b in zip(short[i],short[j])])
    cands.sort(key=lambda v:sum(abs(x) for x in v)); tried=0; nf=0
    for v in cands:
        l1=sum(abs(x) for x in v)
        if l1==0 or l1>maxl1 or tried>=maxcand: continue
        tried+=1; c=v[0]; coef=[(i-1,v[i]) for i in range(1,N) if v[i]]; x=F.norm_fixed(c,coef)
        if x==0 or x%q or (x//q)%q==0: continue
        m=x//q
        while m%F.p==0: m//=F.p
        while True:
            g=math.gcd(m,P)
            if g==1: break
            m//=g
        if m>=q and flint.fmpz(m).is_probable_prime(): continue         # one large prime factor >= q: it is split in almost all cases
        nf+=1; f=pari.factor(x); fac=[(int(f[0][k]),int(f[1][k])) for k in range(len(f[0]))]
        if all(r==q or r==F.p or (not F.split(r)) or r<q for r,e in fac): return (q,c,coef,fac),nf
    return None,nf
if __name__=="__main__":
    p,n,start,cnt=int(sys.argv[1]),int(sys.argv[2]),int(float(sys.argv[3])),int(sys.argv[4]); delta=float(sys.argv[5]) if len(sys.argv)>5 else 0.75
    F,P=setup(p,n); pari=cypari2.Pari(); qs=[]; q=start
    while len(qs)<cnt:
        if F.split(q) and flint.fmpz(q).is_prime(): qs.append(q)
        q+=1
    t=time.time(); ok=0; nf=0
    for q in qs:
        r,k=targeted2(F,P,pari,q,delta); ok+=r is not None; nf+=k
    print("k_{%d^%d}, primes from %.3g, delta=%.2f: %.1f ms per prime, found %d of %d, full factorizations per prime %.2f"%(p,n,start,delta,(time.time()-t)/cnt*1000,ok,cnt,nf/cnt))
