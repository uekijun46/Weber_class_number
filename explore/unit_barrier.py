# Unit (height) method for the cyclotomic Z_2-extension: where can a "too short unit" argument work at all?
# Class r: primes l = 1 mod 2^r exactly; lattice = Z[zeta_{2^r}] * Log(eta_n), rank d = 2^(r-1).
# A hypothetical extra unit eta^(y/l), y in l*L^{-1}, has height |y*Log eta|/l.  Known lower bound for relative units
# (Morisawa-Okazaki): height >= sqrt(N)*log(2+sqrt5), N = 2^n.  Thresholds for l:
#   GH   : typical shortest vector (Gaussian heuristic)  -> below it the extra units are NOT shorter than the bound
#   Blich: worst case with Blichfeldt's bound on the Hermite constant;  Mink: worst case with Minkowski's theorem
import numpy as np, math
from math import lgamma, log, pi, exp
def v0(n):
    N=2**n; M=2**(n+2)
    a=np.empty(N,dtype=np.int64); x=1
    for i in range(N):
        a[i]=x; x=(x*5)%M
    return np.log(np.abs(np.tan(a*math.pi/M)))
L0=log(2+5**.5)
def logV(d): return (d/2)*log(pi)-lgamma(d/2+1)
print(" r   d   n  covol^(1/d)/sqrtN  log10 l_GH  log10 l_Blich  log10 l_Mink")
for r in range(4,10):
    d=2**(r-1)
    for n in (r,r+2,r+5):
        if n>17: continue
        v=v0(n); N=2**n
        V=np.array([np.roll(v,-j*2**(n-r)) for j in range(d)])
        G=V@V.T
        s,ld=np.linalg.slogdet(G)
        lc=ld/(2*d)                      # log covol^(1/d)
        m=log(math.sqrt(N)*L0)
        gh=d*(lc-logV(d)/d-m)
        mink=d*(lc+log(2)-logV(d)/d-m)
        lgam=log(2/pi)+(2/d)*lgamma(2+d/2)       # log of Blichfeldt's bound for the Hermite constant
        bl=d*(lc+lgam/2-m)
        print(f"{r:2d} {d:4d} {n:3d}   {exp(lc)/math.sqrt(N):8.4f}       {gh/log(10):9.1f}   {bl/log(10):9.1f}    {mink/log(10):9.1f}")
