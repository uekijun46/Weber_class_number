import numpy as np, math
from scipy import integrate, special
L=2*10**8
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
P=sieve(L)
def I(c):
    f=lambda x:(1-math.exp(-(x/c)**2)/math.cosh(x/2))/2*(1/math.sinh(x/2)+1/math.cosh(x/2)) if x<600 else 0.0
    a=integrate.quad(f,1e-12,40,limit=400)[0]
    # beyond 40: integrand ~ 2e^{-x/2}
    b=integrate.quad(lambda x:2*math.exp(-x/2),40,np.inf)[0]
    return a+b
def split_primes(p0,n):
    if p0==2:
        M=2**(n+2); r=P%M; return P[(r==1)|(r==M-1)]
    M=p0**(n+1); res=[r for r in range(1,M) if pow(r,p0-1,M)==1]
    return P[np.isin(P%M,res)]
def rd(p0,n):
    if p0==2: return 2**((n+1)-2.0**(-n))
    e=n+1-p0**(-n)-(1-p0**(-n))/(p0-1); return p0**e
C0=math.pi/2+np.euler_gamma+math.log(8*math.pi)
def psum(sp,c,X):
    q=sp[sp<=min(X,L)].astype(float); t=np.log(q)
    s=0.0
    for m in range(1,6):
        s+=np.sum(2*t/q**(m/2)*np.exp(-(m*t/c)**2)/np.cosh(m*t/2))
    return s
def tail(N,c,X):
    if X<=L: return 0.0
    a,b=math.log(L),math.log(X)
    return 4/N*c*math.sqrt(math.pi)/2*(special.erf(b/c)-special.erf(a/c))
def report(p0,n,hmin,Xs,cs=(20,30,40,50,70,100,150,200,300)):
    N=p0**n; r=rd(p0,n); sp=split_primes(p0,n)
    print(f"\nB_({p0},{n}) deg={N} rd={r:.1f} log rd-C0={math.log(r)-C0:.3f} first split primes={sp[:5].tolist()} hmin={hmin}")
    for X in Xs:
        best=None
        for c in cs:
            B=C0-math.log(r)-I(c)+psum(sp,c,X)+tail(N,c,X)
            if B>0:
                hb=2*c*math.sqrt(math.pi)/(N*B)
                if best is None or hb<best[1]: best=(c,hb,B)
        print(f"  X=1e{math.log10(X):.1f}: "+("B<=0 for all c" if best is None else f"c={best[0]} B={best[2]:.3f} h<{best[1]:.1f}"+("  OK" if best[1]<hmin else "")))
if __name__=="__main__":
    sp=split_primes(23,1)
    print("calib B_(23,1): I(50)=%.6f  sum all split<5e7 (c=50)=%.5f  count=%d"%(I(50),psum(sp,50,5*10**7),np.sum(sp<5*10**7)))
    report(23,1,47,[1e7,5e7,1e8])
    report(19,1,191,[1e4,1e5,1e6])
    report(2,6,10**9,[1e5,1e6,1e7,1e8])
    report(29,1,59,[1e8,1e9,1e10,1e11,1e12,1e13])
    report(31,1,311,[1e9,1e10,1e11,1e12,1e13,1e14])
    report(7,2,97,[1e8,1e9,1e10,1e11,1e12,1e13,1e14])
    report(3,4,163,[1e8,1e9,1e10,1e11,1e12,1e14,1e16])
    report(2,7,10**9,[1e12,1e16,1e20,1e24,1e30])
    report(5,3,251,[1e12,1e20,1e30,1e40])
    report(37,1,149,[1e12,1e14,1e16,1e18])
