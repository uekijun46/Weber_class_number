# Exploration (untrusted): how often is the norm of a small element of k_p a product of small split primes?
#   python3 probe.py p count_per_weight
import sys, math, random, time, flint
sys.path.insert(0,"../../relations")
from field import Field
p=int(sys.argv[1]); K=int(sys.argv[2]); F=Field(p,1); F.fixed(320); N=F.N
Ys=[10**4,10**5,10**6]; LP=10**8
sp=[q for q in range(2,Ys[-1]) if F.split(q) and flint.fmpz(q).is_prime()]
P={Y:math.prod(q for q in sp if q<Y) for Y in Ys}
print("k_%d: split primes below 1e4, 1e5, 1e6: %s; smallest %s"%(p,[sum(1 for q in sp if q<Y) for Y in Ys],sp[:6]),flush=True)
rnd=random.Random(7)
for w in (2,3,4,5,6):
    t=time.time(); full={Y:0 for Y in Ys}; part={Y:0 for Y in Ys}; logs=[]; seen=set(); n=0
    tries=0
    while n<K and tries<20*K:
        tries+=1
        idx=tuple(sorted(rnd.sample(range(1,N),w-1))); sg=tuple(rnd.choice((1,-1)) for _ in idx); c=rnd.choice((-2,-1,0,1,2))
        key=(c,idx,sg)
        if key in seen: continue
        seen.add(key); x=F.norm_fixed(c,[(0,1)]+list(zip(idx,sg)))
        if x<2: continue
        n+=1; logs.append(math.log10(x))
        while x%p==0: x//=p
        for Y in Ys:
            m=x
            while True:
                g=math.gcd(m,P[Y])
                if g==1: break
                m//=g
            if m==1: full[Y]+=1
            elif m<LP and flint.fmpz(m).is_prime(): part[Y]+=1
    logs.sort()
    print("weight %d: %d elements, log10|N| median %.1f, min %.1f; all split prime factors below 1e4/1e5/1e6: %s; the same with one extra prime below 1e8: %s  (%.0fs)"%(w,n,logs[len(logs)//2],logs[0],[full[Y] for Y in Ys],[part[Y] for Y in Ys],time.time()-t),flush=True)
