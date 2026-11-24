# Floating point scan (no proof): for which layers k_{p^n} do very small primes split completely, and how much does it help?
import numpy as np, math
C0=math.pi/2+0.5772156649015329+math.log(8*math.pi)
def Gg(c):   # Gaussian F (GRH version)
    x=np.concatenate([np.linspace(1e-9,1,20001),np.linspace(1,60,60001)[1:],np.linspace(60,3000,60001)[1:]])
    F=np.exp(-(x/c)**2); y=(1-F)/2*(1/np.sinh(np.minimum(x/2,300))+1/np.cosh(np.minimum(x/2,300)))
    return float(np.trapezoid(y,x))
LIM=2*10**6; sv=np.ones(LIM+1,bool); sv[:2]=False
for i in range(2,int(LIM**.5)+1):
    if sv[i]: sv[i*i::i]=False
primes=np.nonzero(sv)[0]
def hmin(N,p):
    best=None
    for l in primes:
        l=int(l)
        if best is not None and l>=best: break
        if l==p: continue
        f=1; x=l%N
        while x!=1: x=x*l%N; f+=1
        if best is None or l**f<best: best=l**f
    return best
cs=[5,6,7,8,9,10,12]; Gc={c:Gg(c) for c in cs}
fields=[(p,1) for p in [int(x) for x in primes if 31<=x<200]]+[(3,5),(5,3),(7,3),(11,2),(13,2)]
rows=[]
for p,n in fields:
    N=p**n; M=p**(n+1); ex=sum((p-1)*p**(j-1)*(j+1) for j in range(1,n+1)); logrd=ex/N*math.log(p); hm=hmin(N,p)
    res=np.array([pow(r,p-1,M)==1 for r in range(M)]); sp=[int(q) for q in primes[res[primes%M]]]
    # under GRH: least number of smallest split primes with h < hm
    t_need=None
    for t in range(1,min(len(sp),400)+1):
        ok=False
        for c in cs:
            S=sum(2*math.log(q)/q**(m/2)*math.exp(-(m*math.log(q)/c)**2) for q in sp[:t] for m in range(1,12))
            B=C0-logrd-Gc[c]+S
            if B>0 and 2*c*math.sqrt(math.pi)*math.exp((c/4)**2)/(N*B)<hm: ok=True; break
        if ok: t_need=t; break
    # without GRH: sum up to 2e6 (c large), then the slope 4/N per e-fold
    S0=sum(4*math.log(q)/(q**m+1) for q in sp for m in (1,2,3)); need=logrd+1.2639-C0
    lgX=(math.log(LIM)+max(0,need-S0)*N/4)/math.log(10)
    rows.append((p,n,N,math.exp(logrd),hm,sp[:4],t_need,sp[t_need-1] if t_need else None,lgX))
print("field | degree | rd | threshold | smallest split primes | GRH: primes needed (largest) | no GRH: log10 of the range needed")
for p,n,N,rd,hm,s4,t,qt,lgX in rows:
    print("k_{%d^%d} | %d | %.0f | %d | %s | %s (%s) | %.0f"%(p,n,N,rd,hm,s4,t,qt,lgX))
