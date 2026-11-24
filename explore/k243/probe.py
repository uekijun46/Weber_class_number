# Exploration (untrusted, floating/fixed point): can relations over a factor base settle the smallest split primes of k_{3^5}?
# A prime that does not split completely in k_{3^5} is settled if h(k_{3^4}) = 1, so only the split primes (q = +-1 mod 729)
# matter.  An element is called smooth for the bound Y if every split prime factor of its norm is below Y.
#   python3 probe.py [number of weight-3 elements]
import sys, math, random, time, itertools, flint
sys.path.insert(0,"../../relations")
from field import Field
F=Field(3,5); F.fixed(360); N=F.N; M=F.M
Ys=[10**5,10**6,10**7,10**8]; YM=Ys[-1]
t=time.time(); sv=bytearray([1])*(YM+1); sv[0]=sv[1]=0
for i in range(2,int(YM**.5)+1):
    if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
split=[q for r in (1,M-1) for q in range(r,YM,M) if q>1 and sv[q]]; split.sort()
small=[r for r in range(2,3000) if sv[r] and r!=3 and r%M not in (1,M-1)]
print("split primes of k_{3^5} below",Ys,":",[sum(1 for q in split if q<Y) for Y in Ys],"(%.0fs)"%(time.time()-t),flush=True)
def classify(x):          # returns the largest split prime factor found if the norm is smooth for 10^8, else None
    while x%3==0: x//=3
    big=1
    for q in split:
        if q*q*q>x and x<YM*YM: break
        while x%q==0: x//=q; big=q
    for r in small:        # non-split primes divide the norm with exponent >= 3
        while x%r==0: x//=r
    if x==1: return big
    if x<YM and x%M in (1,M-1) and sv[x]: return x
    if x<YM*YM:
        # at most two prime factors left; smooth only if both are split primes below 10^8 (or non-split, which is rare)
        for q in split:
            if q*q>x: break
            if x%q==0:
                y=x//q
                if y<YM and sv[y] and y%M in (1,M-1): return max(big,q,y)
    return None
def report(name,items):
    t=time.time(); cnt={Y:0 for Y in Ys}; n=0; logs=[]; triv=[0]
    for c,coef in items:
        x=F.norm_fixed(c,coef)
        if x<2: continue
        n+=1; logs.append(math.log10(x)); b=classify(x)
        if b is not None:
            if b==1: triv[0]+=1                       # no split prime at all: such an element gives no information on split primes
            else:
                for Y in Ys:
                    if b<Y: cnt[Y]+=1
    logs.sort()
    print("%s: %d elements, median log10|N| = %.1f; without any split prime factor: %d; with split prime factors, all below Y = 1e5, 1e6, 1e7, 1e8: %s (%.0fs)"%(name,n,logs[len(logs)//2] if logs else 0,triv[0],[cnt[Y] for Y in Ys],time.time()-t),flush=True)
w2=[(c,[(0,1),(j,s)]) for c in (-2,-1,0,1,2) for j in range(1,N) for s in (1,-1)]
report("weight 2 (all)",w2)
rnd=random.Random(1); K3=int(sys.argv[1]) if len(sys.argv)>1 else 20000
w3=[]
for _ in range(K3):
    j,k=sorted(rnd.sample(range(1,N),2)); w3.append((rnd.choice((-2,-1,0,1,2)),[(0,1),(j,rnd.choice((1,-1))),(k,rnd.choice((1,-1)))]))
report("weight 3 (sample)",w3)
# elements lying in a prime above a small split prime q: short vectors of the lattice of index q
tg=[]
for q in split[:0]:
    R=F.periods_mod(q); B=flint.fmpz_mat(N,N); B[0,0]=q
    for i in range(N-1): B[i+1,0]=(-R[i])%q; B[i+1,i+1]=1
    L=B.lll(); vs=[[int(L[i,j]) for j in range(N)] for i in range(N)]
    vs.sort(key=lambda v:sum(abs(x) for x in v))
    for v in vs[:80]:
        if 0<sum(abs(x) for x in v)<=6: tg.append((v[0],[(i-1,v[i]) for i in range(1,N) if v[i]]))
report("short vectors in the primes above the 40 smallest split primes",tg)
