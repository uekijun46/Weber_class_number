# Prototype for k_{7^2}: prove principality of the primes above split primes from norms of small elements only.
# A row (alpha, factorization of |N(alpha)|) certifies a split prime q if q divides the norm exactly once and every other
# prime factor is 7, a non-split prime, or an already certified split prime.
import sys, math, random, time, pickle, os
from multiprocessing import Pool
from mpmath import mp, mpf, cos, pi
import cypari2
P0,NN,M=7,49,343
def periods(S=170):
    g=next(x for x in range(2,M) if x%7 and pow(x,294//2,M)!=1 and pow(x,294//3,M)!=1 and pow(x,294//7,M)!=1)
    H=[pow(g,NN*t,M) for t in range(6)]; mp.prec=240
    eta=[sum(cos(2*pi*(pow(g,i,M)*h%M)/M) for h in H) for i in range(NN)]
    return [int(e*mpf(2)**S) for e in eta], S
E,S=periods()
def norm(c,idx,sg):
    prod=1
    for j in range(NN):
        v=c<<S
        for i,s in zip(idx,sg): v+=s*E[(i+j)%NN]
        prod*=v
    sh=S*NN; prod=abs(prod); return (prod+(1<<(sh-1)))>>sh
def work(args):
    seed,n=args; rnd=random.Random(seed); pari=cypari2.Pari(); out=[]; t0=time.time(); tf=0
    for _ in range(n):
        w=rnd.choice((1,2,2,3,3,3,4,4,5)); idx=sorted(rnd.sample(range(NN),w)); idx=[i-idx[0] for i in idx]
        sg=[rnd.choice((-1,1)) for _ in idx]; c=rnd.choice((-2,-1,0,0,1,2))
        x=norm(c,idx,sg)
        if x<2 or x>10**26: continue
        t=time.time(); f=pari.factor(x); tf+=time.time()-t
        fac=[(int(f[0][k]),int(f[1][k])) for k in range(len(f[0]))]
        out.append((c,tuple(idx),tuple(sg),tuple(fac)))
    return out,time.time()-t0,tf
def split(q): return q!=7 and pow(q,6,M)==1
def peel(rows):
    cert={}; order=[]                         # q -> row index
    need=[]; byprime={}
    for r,(c,idx,sg,fac) in enumerate(rows):
        sp=[(q,e) for q,e in fac if split(q)]
        need.append(dict(sp))
        for q,e in sp: byprime.setdefault(q,[]).append(r)
    pending=[set(d) for d in need]; stack=[r for r in range(len(rows)) if len(pending[r])==1]
    while stack:
        r=stack.pop()
        if len(pending[r])!=1: continue
        (q,)=pending[r]
        if need[r][q]!=1 or q in cert: 
            if q in cert: pending[r].discard(q)
            continue
        cert[q]=r; order.append(q)
        for r2 in byprime[q]:
            if q in pending[r2]:
                pending[r2].discard(q)
                if len(pending[r2])==1: stack.append(r2)
    return cert,order
if __name__=="__main__":
    n=int(sys.argv[1]); seed0=int(sys.argv[2]); X0=int(sys.argv[3])
    rows=pickle.load(open("rows.pkl","rb")) if os.path.exists("rows.pkl") else []
    t0=time.time()
    with Pool(2) as pool: res=pool.map(work,[(seed0,n//2),(seed0+1,n//2)])
    for out,tt,tf in res: rows+=out
    print("elements this run: %d, kept rows: %d total, %.0f s wall; per element %.0f us (of which factoring %.0f us)"%(n,len(rows),time.time()-t0,res[0][1]/(n//2)*1e6,res[0][2]/(n//2)*1e6))
    pickle.dump(rows,open("rows.pkl","wb"))
    t1=time.time(); cert,order=peel(rows)
    sv=bytearray([1])*(X0+1)
    for i in range(2,int(X0**.5)+1):
        if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
    sp=[q for q in range(2,X0+1) if sv[q] and split(q)]
    miss=[q for q in sp if q not in cert]
    print("peeling %.1f s: certified primes: %d (of which > X0: %d); split primes <= %d: %d, missing: %d"%(time.time()-t1,len(cert),sum(1 for q in cert if q>X0),X0,len(sp),len(miss)))
    for b in (10**4,10**5,2*10**5,5*10**5,10**6):
        if b<=X0: print("   <= %d: %d of %d certified"%(b,sum(1 for q in sp if q<=b and q in cert),sum(1 for q in sp if q<=b)))
    print("   smallest missing:",miss[:12])
