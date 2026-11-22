import mpmath as mp
from sympy import factorint, primerange, n_order
mp.mp.dps=30
def logdisc(m):
    # conductor-discriminant formula
    L=mp.mpf(0)
    for p,v in factorint(m).items():
        ex = sum((p**k-p**(k-1))*((k+2) if p==2 else (k+1)) for k in range(1,v+1))
        L += (m//p**v)*ex*mp.log(p)
    return L
def frob_order(q,m):
    # order of Frobenius of q (q not dividing m) in Gal(Q(m)/Q)
    o=1
    for p,v in factorint(m).items():
        if p==2:
            M=2**(v+2); k=0
            while True:
                e=pow(q,2**k,M)
                if e in (1,M-1): break
                k+=1
            o=max(o,2**k) if True else 0
        else:
            M=p**(v+1); k=0
            while pow(q,(p-1)*p**k,M)!=1: k+=1
            o=max(o,p**k)
    # orders are prime powers for distinct primes: lcm = product
    return o
def frob_order_lcm(q,m):
    from math import lcm
    o=1
    for p,v in factorint(m).items():
        if p==2:
            M=2**(v+2); k=0
            while pow(q,2**k,M) not in (1,M-1): k+=1
            o=lcm(o,2**k)
        else:
            M=p**(v+1); k=0
            while pow(q,(p-1)*p**k,M)!=1: k+=1
            o=lcm(o,p**k)
    return o
def principal_primes(m, Q, split_cert=()):
    """list of (norm, count) of prime ideals certified principal via A1/A2, plus certified totally split primes"""
    S=[]
    fm=factorint(m)
    for q in primerange(2,Q+1):
        if q in fm:
            v=fm[q]; mp_=m//q**v
            f=frob_order_lcm(q,mp_) if mp_>1 else 1
            if f>1:
                g=m//(q**v*f); S.append((q**f,g))
        else:
            f=frob_order_lcm(q,m)
            if f>1: S.append((q**f,m//f))
    for q in split_cert: S.append((q,m))
    return S
I_cache={}
def integral(c):
    F=lambda x: mp.e**(-(x/c)**2)/mp.cosh(x/2)
    g=lambda x: (1-F(x))/2*(1/mp.sinh(x/2)+1/mp.cosh(x/2))
    return mp.quad(g,[0,1,5,20,mp.inf])
def Bval(m,c,S):
    N=m
    F=lambda x: mp.e**(-(x/c)**2)/mp.cosh(x/2)
    lrd=logdisc(m)/N
    s=mp.mpf(0)
    for (nP,cnt) in S:
        lp=mp.log(nP); k=1
        while True:
            t=lp/mp.power(nP,mp.mpf(k)/2)*F(k*lp)
            if t<mp.mpf(10)**-25: break
            s+=cnt*t; k+=1
    B=mp.pi/2+mp.euler+mp.log(8*mp.pi)-lrd-integral(c)+2*s/N
    return B
def bound(m,c,S):
    B=Bval(m,c,S)
    return (2*c*mp.sqrt(mp.pi)/(m*B) if B>0 else mp.inf), B
def best(m,S,cs):
    res=[(bound(m,c,S)[0],c) for c in cs]
    return min(res)
def hmin(m):
    best=None
    fm=factorint(m)
    for l in primerange(2,10**4):
        if l in fm: continue
        f=n_order(l,m) if m>1 else 1
        v=l**f
        if best is None or v<best: best=v
    return best
