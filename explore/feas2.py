import numpy as np, math
from scipy import integrate, special
import feas
from feas import *
def Igrh(c):
    f=lambda x:(1-math.exp(-(x/c)**2))/2*(1/math.sinh(x/2)+1/math.cosh(x/2))
    return integrate.quad(f,1e-12,60,limit=400)[0]+integrate.quad(lambda x:2*math.exp(-x/2),60,np.inf)[0]
def ram(p0,N,c,grh=False):
    s=0
    for m in range(1,60):
        t=m*math.log(p0); F=math.exp(-(t/c)**2)/(1 if grh else math.cosh(t/2))
        s+=2/N*math.log(p0)/p0**(m/2)*F
    return s
def psum_grh(plist,c):
    s=0
    for p in plist:
        for m in range(1,40):
            t=m*math.log(p); s+=2*math.log(p)/p**(m/2)*math.exp(-(t/c)**2)
    return s
# validate GRH formula on conductor 512 (Miller: c=8.7, ten primes, contribution>0.6898, h<147)
S=[3583,5119,6143,7681,8191,10753,11777,12289,12799,13313]
c=8.7; r=rd(2,7); ps=psum_grh(S,c); B=C0-math.log(r)-Igrh(c)+ps
print("512 GRH check: primesum=%.4f B=%.4f h<%.1f"%(ps,B,2*c*math.sqrt(math.pi)*math.exp((c/4)**2)/(128*B)))
# B_(7,2) under GRH with few primes
sp=split_primes(7,2); r=rd(7,2); N=49
print("B_(7,2) split primes:",sp[:12].tolist())
for k in (0,1,2,3,5,10):
    best=None
    for c in np.arange(2,14,0.1):
        B=C0-math.log(r)-Igrh(c)+psum_grh(sp[:k].tolist(),c)+ram(7,N,c,True)
        if B>0:
            hb=2*c*math.sqrt(math.pi)*math.exp((c/4)**2)/(N*B)
            if best is None or hb<best[1]: best=(c,hb,B)
    print("  GRH, first %d split primes + ramified 7:"%k, None if best is None else "c=%.1f B=%.3f h<%.1f"%(best[0],best[2],best[1]))
# unconditional refined, with ramified prime, with/without 19
def uncond(p0,n,Xs,drop=(),cs=(20,30,40,50,60,70,85,100,150,200,300),hmin=None):
    N=p0**n; r=rd(p0,n); sp=split_primes(p0,n); sp=sp[~np.isin(sp,drop)]
    for X in Xs:
        best=None
        for c in cs:
            B=C0-math.log(r)-I(c)+psum(sp,c,X)+tail(N,c,X)+ram(p0,N,c)
            if B>0:
                hb=2*c*math.sqrt(math.pi)/(N*B)
                if best is None or hb<best[1]: best=(c,hb,B)
        cnt=np.sum(sp<=X) if X<=feas.L else X/math.log(X)/N
        print(f"  B_({p0},{n}) X={X:.0e} (#split~{cnt:.2e}) drop={drop}: "+("B<=0" if best is None else f"c={best[0]} B={best[2]:.3f} h<{best[1]:.1f}"))
print("unconditional:")
uncond(7,2,[1e6,3e6,1e7,3e7,1e8,2e8])
uncond(7,2,[1e8,1e10],drop=(19,))
uncond(29,1,[1e9,3e9,1e10])
uncond(3,4,[1e10,3e10,1e11,1e12])
uncond(31,1,[1e11,3e11,1e12])
uncond(2,7,[1e22,1e23,1e24])
uncond(11,2,[1e12,1e20,1e30])
uncond(13,2,[1e20,1e40])
