# Untrusted generator: for every split prime q in [A,B) find alpha with q || N(alpha) such that every other split prime factor is < q.
#   python3 gen_range.py p n A B OUT [procs] [time limit]
import sys, time, random
from multiprocessing import Pool
import flint, cypari2
from field import Field
_F=None; _pari=None
def _init(p,n):
    global _F,_pari
    _F=Field(p,n); _F.fixed(); _pari=cypari2.Pari()
def one(q):
    F=_F; N=F.N; p=F.p; rnd=random.Random(q)
    for attempt in range(6):
        R=F.periods_mod(q,rnd)
        B=flint.fmpz_mat(N,N); B[0,0]=q
        perm=list(range(N-1))
        if attempt: rnd.shuffle(perm)                     # a different basis order gives different short vectors
        for k,i in enumerate(perm): B[k+1,0]=(-R[i])%q; B[k+1,k+1]=1
        L=B.lll(); vs=[]
        for i in range(N):
            v=[0]*N; v[0]=int(L[i,0])
            for k,j in enumerate(perm): v[j+1]=int(L[i,k+1])
            vs.append(v)
        vs.sort(key=lambda v:sum(abs(x) for x in v)); short=vs[:16+8*attempt]; cands=list(vs)
        for i in range(len(short)):
            for j in range(i+1,len(short)):
                cands.append([a+b for a,b in zip(short[i],short[j])]); cands.append([a-b for a,b in zip(short[i],short[j])])
        cands.sort(key=lambda v:sum(abs(x) for x in v))
        for v in cands[:600+400*attempt]:
            l1=sum(abs(x) for x in v)
            if l1==0 or l1>20: continue
            c=v[0]; coef=[(i-1,v[i]) for i in range(1,N) if v[i]]
            x=F.norm_fixed(c,coef)
            if x==0 or x%q or (x//q)%q==0: continue
            f=_pari.factor(x//q); fac=[(int(f[0][k]),int(f[1][k])) for k in range(len(f[0]))]
            if all(r==p or (not F.split(r)) or r<q for r,e in fac):
                fac=sorted(fac+[(q,1)])
                return "%d | %d | %s | %s"%(q,c," ".join("%d:%d"%ia for ia in coef)," ".join("%d^%d"%re for re in fac))
    return "FAIL %d"%q
if __name__=="__main__":
    p,n,A,B=[int(x) for x in sys.argv[1:5]]; outfn=sys.argv[5]; procs=int(sys.argv[6]) if len(sys.argv)>6 else 2
    F=Field(p,n); t0=time.time(); M=F.M              # no PARI in the parent process (it must not be initialized before the fork)
    res=sorted(r for r in range(M) if r%p and pow(r,p-1,M)==1); qs=[]
    base=A-(A%M)
    while base<B:
        for r in res:
            q=base+r
            if A<=q<B and flint.fmpz(q).is_probable_prime(): qs.append(q)
        base+=M
    out=open(outfn,"w"); out.write("# field %d %d\n"%(p,n)); fails=[]
    with Pool(procs,initializer=_init,initargs=(p,n)) as pool:
        for line in pool.imap(one,qs,chunksize=32):
            if line.startswith("FAIL"): fails.append(int(line.split()[1]))
            else: out.write(line+"\n")
    out.close()
    print("range [%d,%d): %d split primes, %d rows, %d failures %s, %.0fs (%.0f ms per prime per core)"%(A,B,len(qs),len(qs)-len(fails),len(fails),fails[:5],time.time()-t0,(time.time()-t0)*procs/max(1,len(qs))*1000))
