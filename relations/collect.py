# Random small elements and the factorizations of their norms (untrusted search step).
#   python3 collect.py p n count seed [maxweight] [normlimit_digits]
import sys, os, time, pickle, random
from multiprocessing import Pool
import cypari2
from field import Field
def work(args):
    p,n,seed,cnt,maxw,lim=args; F=Field(p,n); F.fixed(); rnd=random.Random(seed); pari=cypari2.Pari(); out=[]; N=F.N
    for _ in range(cnt):
        w=rnd.randint(1,maxw); idx=sorted(rnd.sample(range(N),w)); idx=[i-idx[0] for i in idx]
        coef=[(i,rnd.choice((-1,1))) for i in idx]; c=rnd.choice((-2,-1,0,0,1,2))
        x=F.norm_fixed(c,coef)
        if x<2 or x>10**lim: continue
        f=pari.factor(x); out.append((c,coef,tuple((int(f[0][k]),int(f[1][k])) for k in range(len(f[0])))))
    return out
if __name__=="__main__":
    p,n,cnt,seed=[int(x) for x in sys.argv[1:5]]; maxw=int(sys.argv[5]) if len(sys.argv)>5 else 5; lim=int(sys.argv[6]) if len(sys.argv)>6 else 26
    fn="rows_%d_%d.pkl"%(p,n); rows=pickle.load(open(fn,"rb")) if os.path.exists(fn) else []; t0=time.time()
    with Pool(2) as pool: res=pool.map(work,[(p,n,seed,cnt//2,maxw,lim),(p,n,seed+1,cnt//2,maxw,lim)])
    for o in res: rows+=o
    pickle.dump(rows,open(fn,"wb")); print("%d elements in %.0fs; rows kept in total: %d"%(cnt,time.time()-t0,len(rows)))
