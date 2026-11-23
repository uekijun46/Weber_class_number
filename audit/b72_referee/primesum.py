# lower bound for 2 sum_{p in S} sum_{m>=1} log p p^{-m/2} F(m log p), using mpmath interval arithmetic
import sys
from mpmath import iv, mp
iv.prec=100
P=[int(l) for l in open(sys.argv[1])]
cs=[int(a) for a in sys.argv[2:]]
for c in cs:
    C=iv.mpf(c); S=iv.mpf(0)
    for p in P:
        L=iv.log(p)
        for m in range(1,11):
            x=m*L
            S+= 2*L*iv.exp(-(x/C)**2)/(iv.mpf(p)**m+1)   # p^{-m/2}/cosh(m log p/2)=2/(p^m+1)
    S=2*S
    print(c, mp.nstr(S.a,20), mp.nstr(S.b,20), flush=True)
