# How many split primes does the unconditional method need for k_{3^4}?  (floating point estimate)
import numpy as np, math, mpmath as mp, time
M=243; N=81; logrd=364/81*math.log(3); C0=math.pi/2+0.5772156649015329+math.log(8*math.pi); hmin=163
cs=[60,80,100,150,200,300]
def I(c):
    F=lambda x: mp.exp(-(x/c)**2)/mp.cosh(x/2)
    return float(mp.quad(lambda x:(1-F(x))/2*(1/mp.sinh(x/2)+1/mp.cosh(x/2)),[0,1,10,50,200,2000]))
Ic={c:I(c) for c in cs}
XMAX=10**11; SEG=4*10**7                      # segment length in k (q = 243k +- 1)
small=[r for r in range(2,int(XMAX**.5)+2) if all(r%d for d in range(2,int(r**.5)+1))] if False else None
lim=int(XMAX**.5)+1; sv=np.ones(lim+1,bool); sv[:2]=False
for i in range(2,int(lim**.5)+1):
    if sv[i]: sv[i*i::i]=False
small=np.nonzero(sv)[0]; small=small[small!=3]
S={c:0.0 for c in cs}; cnt=0; t0=time.time(); checkpoints=[10**8,10**9,2*10**9,5*10**9,10**10,2*10**10,3*10**10,5*10**10,10**11]; ci=0
kmax=XMAX//M+1; k0=1
inv=np.array([pow(243,-1,int(r)) for r in small],dtype=np.int64)
while k0<=kmax and ci<len(checkpoints):
    k1=min(k0+SEG,checkpoints[ci]//M+1); L=k1-k0
    for sign in (1,-1):
        a=np.ones(L,bool)
        qhi=243*(k1-1)+sign
        for r,iv in zip(small,inv):
            r=int(r)
            if r*r>qhi: break
            kk=(-sign*iv)%r                      # 243 k + sign = 0 mod r
            st=(kk-k0)%r
            if 243*(k0+st)+sign==r: st+=r
            a[st::r]=False
        q=243*(k0+np.nonzero(a)[0]).astype(np.float64)+sign
        q=q[q>1]; cnt+=len(q); lq=np.log(q)
        for c in cs:
            S[c]+=np.sum(4*lq*np.exp(-(lq/c)**2)/(q+1))+np.sum(2*lq/q*np.exp(-(2*lq/c)**2)/np.cosh(lq))
    k0=k1
    if 243*k0>=checkpoints[ci]:
        X=checkpoints[ci]; ci+=1; best=None
        for c in cs:
            B=C0-logrd-Ic[c]+S[c]
            if B>0:
                h=2*c*math.sqrt(math.pi)/(N*B)
                if best is None or h<best[0]: best=(h,c,B)
        print("X=%.0e: %d split primes; "%(X,cnt)+("best h < %.1f (c=%d, B=%.4f)"%best if best else "B <= 0 for all c; max B = %.4f"%max(C0-logrd-Ic[c]+S[c] for c in cs))+"  [%.0fs]"%(time.time()-t0),flush=True)
        if best and best[0]<hmin*0.6: break
