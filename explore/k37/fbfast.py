# Exploration (untrusted search): relations for the first layer k_p from small elements, fast version.
# Elements: c + sum_{i in I} s_i eta_i with |I| = w, s_i = +-1, c in -2..2, one element for each orbit under Gal(k_p/Q) and sign.
# The norm is computed modulo six primes l = 1 mod p^2 with numpy; an element is kept if its norm is
#   p^a * (split primes below Y) * (at most one more prime below LP).
#   python3 fbfast.py p w J:j Y LP OUT      (the part j of J of the supports)
import sys, math, itertools, pickle, time, numpy as np, flint
sys.path.insert(0,"../../relations")
from field import Field
p,w=int(sys.argv[1]),int(sys.argv[2]); J,j0=[int(u) for u in sys.argv[3].split(":")]; Y,LP,out=int(float(sys.argv[4])),int(float(sys.argv[5])),sys.argv[6]
F=Field(p,1); N=F.N; M=F.M
ls=[]; l=(2**31//M)*M+1
while len(ls)<6:
    l-=M
    if flint.fmpz(l).is_prime(): ls.append(l)
Rl=[np.array(F.periods_mod(l),dtype=np.int64) for l in ls]
Mc=math.prod(ls); cf=[(Mc//l)*pow(Mc//l,-1,l) for l in ls]
fb=[q for q in range(2,Y) if F.split(q) and flint.fmpz(q).is_prime()]; P=flint.fmpz(math.prod(fb)); Pi=math.prod(fb)
sg=np.array([(1,)+s for s in itertools.product((1,-1),repeat=w-1)],dtype=np.int64); cs=np.array([-2,-1,0,1,2],dtype=np.int64)
def canonical(I):                          # I: sorted, contains 0; True if I is the smallest among its translates that contain 0
    for a in I[1:]:
        if tuple(sorted((i-a)%N for i in I))<I: return False
    return True
def supports():
    for t,idx in enumerate(itertools.combinations(range(1,N),w-1)):
        if t%J==j0 and canonical((0,)+idx): yield (0,)+idx
rels=[]; n=0; t0=time.time(); tn=tp=0
def batch(Is):
    global n,tn,tp
    t=time.time(); I=np.array(Is,dtype=np.int64); res=[]
    for l,R in zip(ls,Rl):
        prod=np.ones((len(Is),len(sg),len(cs)),dtype=np.int64)
        for j in range(N):
            v=R[(I+j)%N]@sg.T
            prod=prod*((v[:,:,None]+cs[None,None,:])%l)%l
        res.append(prod.ravel())
    acc=res[0].astype(object)*cf[0]
    for r,c in zip(res[1:],cf[1:]): acc=acc+r.astype(object)*c
    acc=(acc%Mc).tolist(); tn+=time.time()-t; t=time.time(); n+=len(acc); half=Mc//2
    acc=[x if x<=half else Mc-x for x in acc]
    for k0 in range(0,len(acc),256):
        blk=acc[k0:k0+256]; lev=[[flint.fmpz(x) if x>1 else flint.fmpz(3) for x in blk]]
        while len(lev[-1])>1:
            a=lev[-1]; lev.append([a[i]*a[i+1] if i+1<len(a) else a[i] for i in range(0,len(a),2)])
        r=[P%lev[-1][0]]
        for a in reversed(lev[:-1]): r=[r[i//2]%a[i] for i in range(len(a))]
        for k,(x,rr) in enumerate(zip(blk,r),k0):
            if x<2: continue
            g=math.gcd(x,int(rr))
            if g==1: continue
            m=x
            while m%p==0: m//=p
            while g>1: m//=g; g=math.gcd(m,g)
            if m!=1 and not (m<LP and flint.fmpz(m).is_prime()): continue
            b,rem=divmod(k,len(sg)*len(cs)); si,c=divmod(rem,len(cs))
            rels.append((int(cs[c]),[(int(i),int(a)) for i,a in zip(Is[b],sg[si])],x,m if m!=1 else None))
    tp+=time.time()-t
B=max(1,200000//(len(sg)*len(cs))); cur=[]
for I in supports():
    cur.append(I)
    if len(cur)==B: batch(cur); cur=[]
if cur: batch(cur)
pickle.dump(rels,open(out,"wb"))
print("k_%d, weight %d, part %d of %d: %d elements in %.0fs (norms %.0fs, tests %.0fs); kept %d (full %d, one large prime %d)"%(p,w,j0,J,n,time.time()-t0,tn,tp,len(rels),sum(1 for r in rels if r[3] is None),sum(1 for r in rels if r[3] is not None)),flush=True)
