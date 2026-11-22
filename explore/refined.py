# Refined bound.  zeta_H = zeta_K * (entire), so the zeros of zeta_K are among those of zeta_H.  Subtracting the two explicit formulas:
#     (h-1) * N * B_cert  <=  W_uncert ,     W_uncert = total weight of the primes of K not known to be principal.
# Non-split primes of K are extended from a lower layer (class number 1), hence principal; so W_uncert only involves split primes > X.
# Compare with Miller's lemma:  h * N * B_cert < 2c sqrt(pi).
import numpy as np, math
from scipy import special
import feas
from feas import split_primes, rd, C0, I, psum, tail
def ram(p0,N,c):
    return sum(2/N*math.log(p0)/p0**(m/2)*math.exp(-(m*math.log(p0)/c)**2)/math.cosh(m*math.log(p0)/2) for m in range(1,60))
SAFETY=2.5       # Brun-Titchmarsh-type allowance for a rigorous upper bound of the tail
def bounds(p0,n,X,sp,cs):
    N=p0**n; r=rd(p0,n); bm=math.inf; br=math.inf
    for c in cs:
        B=C0-math.log(r)-I(c)+psum(sp,c,X)+tail(N,c,X)+ram(p0,N,c)
        if B<=0: continue
        bm=min(bm,2*c*math.sqrt(math.pi)/(N*B))                                   # Miller: h < bm
        W=SAFETY*2*c*math.sqrt(math.pi)*special.erfc(math.log(X)/c)
        br=min(br,W/(N*B))                                                        # refined: h-1 <= br
    return bm,br
def minX(p0,n,sp,cs,f,lo=4.0,hi=13.0):
    if not f(hi): return None
    for _ in range(18):
        mid=(lo+hi)/2
        if f(mid): hi=mid
        else: lo=mid
    return hi
if __name__=="__main__":
    cs=[4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,50,60,80,100,150,200,300]
    print("field      hmin |  Miller h<hmin        | refined h-1<hmin-1    | refined h-1<1 (no rank argument)")
    for (p0,n,hm) in [(7,2,197),(29,1,59),(3,4,163),(31,1,311)]:
        sp=split_primes(p0,n); N=p0**n
        def cnt(x): return np.sum(sp<=10**x) if 10**x<=feas.L else 10**x/(x*math.log(10))/N
        a=minX(p0,n,sp,cs,lambda x: bounds(p0,n,10**x,sp,cs)[0]<hm)
        b=minX(p0,n,sp,cs,lambda x: bounds(p0,n,10**x,sp,cs)[1]<hm-1)
        c=minX(p0,n,sp,cs,lambda x: bounds(p0,n,10**x,sp,cs)[1]<1)
        print(f"B_({p0},{n})  {hm:4d} |  X=1e{a:.2f} (#{cnt(a):.2e}) | X=1e{b:.2f} (#{cnt(b):.2e}) | X=1e{c:.2f} (#{cnt(c):.2e})")
