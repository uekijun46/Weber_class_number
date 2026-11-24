# floating point estimate only: split primes needed for k_p (first layer)
import numpy as np, math, mpmath as mp, sys
p=int(sys.argv[1]); XMAX=int(float(sys.argv[2])); hmin=int(sys.argv[3])
M=p*p; N=p; logrd=(2*(p-1))/p*math.log(p); C0=math.pi/2+0.5772156649015329+math.log(8*math.pi)
cs=[20,30,40,60,80,100,150]
def I(c):
    F=lambda x: mp.exp(-(x/c)**2)/mp.cosh(x/2)
    return float(mp.quad(lambda x:(1-F(x))/2*(1/mp.sinh(x/2)+1/mp.cosh(x/2)),[0,1,10,50,200,2000]))
Ic={c:I(c) for c in cs}
tab=np.zeros(M,bool)
for r in range(1,M):
    if r%p and pow(r,p-1,M)==1: tab[r]=True
lim=int(XMAX**.5)+1; sv=np.ones(lim+1,bool); sv[:2]=False
for i in range(2,int(lim**.5)+1):
    if sv[i]: sv[i*i::i]=False
small=np.nonzero(sv)[0]
S={c:0.0 for c in cs}; cnt=0; first=[]
cps=[x for x in [1e6,1e7,3e7,1e8,3e8,1e9,2e9,4e9,1e10] if x<=XMAX]; ci=0; SEG=5*10**7; lo=2
print("k_%d: log rd = %.4f, C0 - log rd = %.4f, hmin = %d"%(p,logrd,C0-logrd,hmin))
while lo<XMAX and ci<len(cps):
    hi=int(min(lo+SEG,cps[ci])); a=np.ones(hi-lo,bool)
    for r in small:
        r=int(r)
        if r*r>=hi: break
        st=max(r*r,((lo+r-1)//r)*r)-lo
        a[st::r]=False
    q=lo+np.nonzero(a)[0]; q=q[tab[q%M]]
    if len(first)<5: first+=list(q[:5-len(first)])
    cnt+=len(q); qf=q.astype(np.float64); L=np.log(qf)
    for c in cs:
        for m in (1,2,3): S[c]+=float(np.sum(4*L*np.exp(-(m*L/c)**2)/(qf**m+1)))
    lo=hi
    if lo>=cps[ci]:
        best=min(((2*c*math.sqrt(math.pi)/(N*(C0-logrd-Ic[c]+S[c])),c) for c in cs if C0-logrd-Ic[c]+S[c]>0),default=None)
        print("X=%.0e: %d split primes (first %s), "%(cps[ci],cnt,first)+("B<=0" if best is None else "h < %.1f (c=%d)"%best)+", S(c=60)=%.4f"%S[60],flush=True); ci+=1
