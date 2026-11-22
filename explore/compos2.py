# Composita K = prod B_{p_i,n_i}: cyclic over Q of degree N.  Assume every proper subfield has class number 1 (induction).
# Free primes (automatically principal): unramified q that is not totally split (extended from the decomposition field);
# ramified p_i unless p_i splits completely in the other factors.  Only the remaining primes need certificates.
import math, itertools, numpy as np
from feas import C0, I
L=2*10**7
sv=np.ones(L+1,bool); sv[:2]=False
for i in range(2,int(L**.5)+1):
    if sv[i]: sv[i*i::i]=False
PR=np.nonzero(sv)[0]
def nchar(p,a): return 1 if a==0 else p**a-p**(a-1)
def cexp(p,a): return 0 if a==0 else (a+2 if p==2 else a+1)
def logrd(F):
    N=math.prod(p**n for p,n in F); s=0
    for a in itertools.product(*[range(n+1) for p,n in F]):
        s+=math.prod(nchar(p,ai) for (p,n),ai in zip(F,a))*sum(cexp(p,ai)*math.log(p) for (p,n),ai in zip(F,a))
    return s/N
def gvec(p,n,qs):       # number of primes of B_{p,n} above each q (q != p), vectorized with python ints via object pow
    M=2**(n+2) if p==2 else p**(n+1)
    r=qs%M; g=np.zeros(len(qs),np.int64)
    # residue classes: compute once per class
    tab={}
    for a in np.unique(r):
        a=int(a)
        if p==2:
            k=0; x=a
            while x%M not in (1,M-1): x=x*x%M; k+=1
        else:
            x=pow(a,p-1,M); k=0
            while x!=1: x=pow(x,p,M); k+=1
        tab[a]=p**(n-k)
    return np.array([tab[int(a)] for a in r],np.int64)
def order(l,N):
    k=1; x=l%N
    while x!=1: x=x*l%N; k+=1
    return k
def setup(F):
    N=math.prod(p**n for p,n in F)
    ram=[p for p,n in F]; mask=~np.isin(PR,ram); qs=PR[mask]
    g=np.ones(len(qs),np.int64)
    for p,n in F: g*=gvec(p,n,qs)
    f=N//g
    rows=[]      # ramified primes
    wief=[]
    for p,n in F:
        gg=1
        for p2,n2 in F:
            if p2!=p: gg*=int(gvec(p2,n2,np.array([p]))[0])
        ff=N//p**n//gg; rows.append((p,gg,ff))
        if gg>1: wief.append((p,gg))
    hm=10**9
    for l in PR[:300]:
        l=int(l)
        if N%l: hm=min(hm,l**order(l,N))
    for p,gg in wief: hm=min(hm,p)
    return N,qs,g,f,rows,hm,wief
def B(F,c,X,S):
    N,qs,g,f,rows,hm,wief=S
    lr=logrd(F); b=C0-lr-I(c)
    use=(g<N)|(qs<=X)
    q=qs[use].astype(float); gg=g[use]; ff=f[use]
    lq=ff*np.log(q)
    sel=lq<8*c; lq=lq[sel]; gg=gg[sel]
    for m in range(1,40):
        t=m*lq; k=t<8*c
        if not k.any(): break
        t=t[k]
        b+=np.sum(2*gg[k]/N*lq[k]*np.exp(-t/2-(t/c)**2)/np.cosh(t/2))
    for p,g1,f1 in rows:
        if f1==1 and p>X: continue
        lq1=f1*math.log(p)
        for m in range(1,200):
            t=m*lq1
            if t>8*c: break
            b+=2*g1/N*lq1*math.exp(-t/2-(t/c)**2)/math.cosh(t/2)
    return b
cs=(6,8,10,12,15,20,25,30,40,50,60,80,100)
def hb(F,X,S):
    N=S[0]; best=math.inf
    for c in cs:
        b=B(F,c,X,S)
        if b>0: best=min(best,2*c*math.sqrt(math.pi)/(N*b))
    return best
if __name__=="__main__":
    Fs=[[(2,1),(3,1)],[(2,2),(3,1)],[(2,1),(5,1)],[(3,1),(5,1)],[(2,1),(7,1)],[(2,1),(3,2)],[(2,2),(5,1)],[(3,1),(7,1)],[(2,3),(3,1)],
        [(2,1),(11,1)],[(2,1),(13,1)],[(2,2),(7,1)],[(2,1),(3,1),(5,1)],[(3,1),(11,1)],[(5,1),(7,1)],[(2,2),(3,2)],[(2,3),(5,1)],[(2,1),(3,1),(7,1)],
        [(3,2),(5,1)],[(2,4),(3,1)],[(2,1),(3,3)],[(2,1),(17,1)],[(3,1),(13,1)],[(2,3),(3,2)],[(2,1),(5,2)],[(2,2),(11,1)],[(2,1),(19,1)]]
    print("field | N | rd | hmin | Wieferich-type (p, #primes) | h-bound with certificates up to X = 0, 1e3, 1e4, 1e5, 1e6, 2e7")
    for F in Fs:
        S=setup(F); N=S[0]
        nsp=[int(np.sum((S[2]==N)&(S[1]<=X))) for X in (1e3,1e4,1e5,1e6,2e7)]
        print(F,N,"%.0f"%math.exp(logrd(F)),S[5],S[6],"|"," ".join(("%.1f"%h if h<1e8 else "-") for h in [hb(F,X,S) for X in (0,1e3,1e4,1e5,1e6,2e7)]),"| #split",nsp,flush=True)
