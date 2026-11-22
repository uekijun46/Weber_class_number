import sys, time
from comp_run import *
from sympy import primerange
def nprimes(Fo,p):
    g=1
    for q,m in Fo: g*=int(compos2.gvec(q,m,np.array([p]))[0])
    return g
res=[]
for Fo in ([(2,1)],[(2,2)],[(3,1)],[(2,1),(3,1)],[(5,1)],[(7,1)],[(3,2)],[(11,1)],[(13,1)]):
    line=[]
    for p in primerange(3,400):
        if p in [q for q,m in Fo]: continue
        t=nprimes(Fo,p)
        if t<2: continue
        ok,rk,t2=genus_ok(sorted(Fo+[(p,1)]),p)
        line.append("%d%s"%(p,"" if ok else "(*t=%d,rk=%d)"%(t2,rk)))
    print("K'=",Fo,": p with >=2 primes:"," ".join(line),flush=True)
