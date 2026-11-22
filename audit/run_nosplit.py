from miller import *
import sys
mp.mp.dps=20
for m in [6,12,10,18,15,14,24,21,20,22,28,36,26]:
    S=principal_primes(m, 20000)
    res=[]
    for c in [2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,28]:
        b,B=bound(m,c,S); res.append((float(b),c,float(B)))
    print(m, hmin(m), min(res), flush=True)
