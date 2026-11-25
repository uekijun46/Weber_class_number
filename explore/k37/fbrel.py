# Exploration (untrusted search): relations over a factor base for the first layer k_p.
# An element c + sum a_i eta_i is kept if its norm is p^a * (split primes below Y) * (at most one more prime below LP).
#   python3 fbrel.py p Y LP weight count seed OUT
import sys, math, random, pickle, time, flint
sys.path.insert(0,"../../relations")
from field import Field
p,Y,LP,w,K,seedarg,out=int(sys.argv[1]),int(float(sys.argv[2])),int(float(sys.argv[3])),int(sys.argv[4]),int(float(sys.argv[5])),sys.argv[6],sys.argv[7]
F=Field(p,1); F.fixed(320); N=F.N
fb=[q for q in range(2,Y) if F.split(q) and flint.fmpz(q).is_prime()]; P=math.prod(fb)
# seed = "J:j" means: all elements of this weight (first index 0, signs +-1, c in -2..2), the part j of J; otherwise a random sample
rels=[]; n=0; t0=time.time()
def gen():
    if ":" in str(seedarg):
        import itertools
        J,j=[int(u) for u in seedarg.split(":")]
        for t,idx in enumerate(itertools.combinations(range(1,N),w-1)):
            if t%J!=j: continue
            for sg in itertools.product((1,-1),repeat=w-1):
                for c in (-2,-1,0,1,2): yield c,idx,sg
    else:
        rnd=random.Random(int(seedarg)); seen=set(); tries=0; got=0
        while got<K and tries<5*K:
            tries+=1
            idx=tuple(sorted(rnd.sample(range(1,N),w-1))); sg=tuple(rnd.choice((1,-1)) for _ in idx); c=rnd.choice((-2,-1,0,1,2))
            if (c,idx,sg) in seen: continue
            seen.add((c,idx,sg)); got+=1; yield c,idx,sg
for c,idx,sg in gen():
    coef=[(0,1)]+list(zip(idx,sg)); x=F.norm_fixed(c,coef)
    if x<2: continue
    n+=1; m=x
    while m%p==0: m//=p
    g=math.gcd(m,P)
    if g==1: continue
    gs=[]
    while g>1: gs.append(g); m//=g; g=math.gcd(m,P)
    if m!=1 and not (m<LP and flint.fmpz(m).is_prime()): continue
    e={}; y=x
    for q in fb:
        if gs[0]%q==0:
            while y%q==0: y//=q; e[q]=e.get(q,0)+1
    rels.append((c,coef,e,m if m!=1 else None))
pickle.dump(rels,open(out,"wb"))
print("k_%d, weight %d: %d elements in %.0fs; kept %d (full %d, with one large prime %d)"%(p,w,n,time.time()-t0,len(rels),sum(1 for r in rels if r[3] is None),sum(1 for r in rels if r[3] is not None)),flush=True)
