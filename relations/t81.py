import sys, time, random, flint, cypari2
from field import Field
F=Field(3,4); F.fixed(); N=F.N; pari=cypari2.Pari()
def targeted(q,maxcand=600,maxl1=18,stats=None):
    t=time.time(); R=F.periods_mod(q); t1=time.time()
    B=flint.fmpz_mat(N,N); B[0,0]=q
    for i in range(N-1): B[i+1,0]=(-R[i])%q; B[i+1,i+1]=1
    L=B.lll(); t2=time.time()
    vs=[[int(L[i,j]) for j in range(N)] for i in range(N)]
    vs.sort(key=lambda v:sum(abs(x) for x in v)); short=vs[:16]; cands=list(vs)
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
        f=pari.factor(x//q); fac=[(int(f[0][k]),int(f[1][k])) for k in range(len(f[0]))]
        if all(r==3 or (not F.split(r)) or r<q for r,e in fac):
            if stats is not None: stats.append((t1-t,t2-t1,time.time()-t2,tried,l1,len(str(x))))
            return (q,c,coef,fac)
    if stats is not None: stats.append((t1-t,t2-t1,time.time()-t2,tried,None,None))
    return None
for start in (10**7,10**9,25*10**9):
    q=start; st=[]; ok=0; n=0; t0=time.time()
    while n<40:
        q=int(pari.nextprime(q+1))
        if not F.split(q): continue
        n+=1; ok+= targeted(q,stats=st) is not None
    import statistics as S
    print("q ~ %.1e: %d/%d found; per prime %.0f ms (periods %.1f, LLL %.0f, candidates %.0f); candidates tried: median %d, max %d; l1 of the element: %s; digits of the norm: %s"%(start,ok,n,(time.time()-t0)/n*1000,S.mean(s[0] for s in st)*1000,S.mean(s[1] for s in st)*1000,S.mean(s[2] for s in st)*1000,S.median(s[3] for s in st),max(s[3] for s in st),sorted(set(s[4] for s in st if s[4])),(min(s[5] for s in st if s[5]),max(s[5] for s in st if s[5]))),flush=True)
