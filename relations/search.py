# Search (untrusted) for rows proving that primes of k_{p^n} above split primes are principal.
# Row: (q, c, coef, fac) with alpha = c + sum a_i eta_i, |N(alpha)| = prod r^e, q || N(alpha), and every other prime factor r
# is p, a non-split prime, or a split prime certified by an earlier row.
import sys, os, time, pickle, random
from multiprocessing import Pool
import flint, cypari2
from field import Field
def peel(F,rows,seeds=()):
    cert={q:None for q in seeds}; order=[]; need=[]; byprime={}
    for r,(c,coef,fac) in enumerate(rows):
        d={q:e for q,e in fac if F.split(q)}; need.append(d)
        for q in d: byprime.setdefault(q,[]).append(r)
    pending=[set(q for q in d if q not in cert) for d in need]
    stack=[r for r in range(len(rows)) if len(pending[r])==1]
    while stack:
        r=stack.pop()
        if len(pending[r])!=1: continue
        (q,)=pending[r]
        if q in cert: pending[r].discard(q); continue
        if need[r][q]!=1: continue
        cert[q]=r; order.append(q)
        for r2 in byprime[q]:
            if q in pending[r2]:
                pending[r2].discard(q)
                if len(pending[r2])==1: stack.append(r2)
    return cert,order
_F=None; _cert=None; _pari=None
def _init(p,n,certset):
    global _F,_cert,_pari
    _F=Field(p,n); _F.fixed(); _cert=certset; _pari=cypari2.Pari()
def targeted(q,maxcand=400,maxl1=16):
    F=_F; N=F.N; R=F.periods_mod(q)
    B=flint.fmpz_mat(N,N)
    B[0,0]=q
    for i in range(N-1): B[i+1,0]=(-R[i])%q; B[i+1,i+1]=1
    L=B.lll(); vs=[[int(L[i,j]) for j in range(N)] for i in range(N)]
    vs.sort(key=lambda v:sum(abs(x) for x in v)); short=vs[:14]; cands=list(vs)
    for i in range(len(short)):
        for j in range(i+1,len(short)):
            cands.append([a+b for a,b in zip(short[i],short[j])]); cands.append([a-b for a,b in zip(short[i],short[j])])
    cands.sort(key=lambda v:sum(abs(x) for x in v)); tried=0
    for v in cands:
        l1=sum(abs(x) for x in v)
        if l1==0 or l1>maxl1 or tried>=maxcand: continue
        tried+=1; c=v[0]; coef=[(i-1,v[i]) for i in range(1,N) if v[i]]
        x=F.norm_fixed(c,coef)
        if x==0 or x%q or (x//q)%q==0: continue
        f=_pari.factor(x); fac=[(int(f[0][k]),int(f[1][k])) for k in range(len(f[0]))]
        if all(r==q or r==F.p or (not F.split(r)) or r in _cert for r,e in fac): return (q,c,coef,fac)
    return (q,None,None,tried)
if __name__=="__main__":
    p,n,X,limit=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]); seeds=[int(x) for x in sys.argv[5].split(",") if x and int(x)] if len(sys.argv)>5 else []
    F=Field(p,n); t0=time.time()
    rows=pickle.load(open("rows_%d_%d.pkl"%(p,n),"rb"))
    extra=pickle.load(open("targeted_%d_%d.pkl"%(p,n),"rb")) if os.path.exists("targeted_%d_%d.pkl"%(p,n)) else []
    cert,order=peel(F,rows+[(c,coef,fac) for (q,c,coef,fac) in extra],seeds)
    sv=bytearray([1])*(X+1)
    for i in range(2,int(X**.5)+1):
        if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
    sp=[q for q in range(2,X+1) if sv[q] and F.split(q)]; miss=[q for q in sp if q not in cert]
    print("split primes <= %d: %d, certified: %d, missing: %d (setup %.0fs)"%(X,len(sp),len(sp)-len(miss),len(miss),time.time()-t0),flush=True)
    maxn=int(sys.argv[6]) if len(sys.argv)>6 else len(miss); miss=miss[:maxn]; fails=[]; done=0; t1=time.time()
    with Pool(2,initializer=_init,initargs=(p,n,set(cert))) as pool:
        for res in pool.imap_unordered(targeted,miss,chunksize=8):
            if res[1] is None: fails.append(res[0])
            else: extra.append(res); done+=1
            if time.time()-t0>limit: pool.terminate(); break
    pickle.dump(extra,open("targeted_%d_%d.pkl"%(p,n),"wb"))
    print("targeted: %d new rows, %d failures in %.0fs (%.0f ms per prime per core); first failures: %s"%(done,len(fails),time.time()-t1,(time.time()-t1)*2/max(1,done+len(fails))*1000,fails[:8]))
