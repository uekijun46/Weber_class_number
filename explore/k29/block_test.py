# Exploration: do small elements of k_p have norms supported only on p and the very small split primes?
import sys, math, itertools
from field import Field
def run(p,n,smalls,wmax=3):
    F=Field(p,n); F.fixed(400); N=F.N; pure=[]; logs=[]; cnt=0
    def test(c,coef):
        nonlocal cnt
        x=F.norm_fixed(c,coef)
        if x<2: return
        cnt+=1; logs.append(math.log10(x)); y=x; e={}
        while y%p==0: y//=p
        for q in smalls:
            while y%q==0: y//=q; e[q]=e.get(q,0)+1
        if y==1 and e: pure.append((math.log10(x),c,coef,e))
    for c in (-2,-1,0,1,2):
        test(c,[(0,1)])
        for j in range(1,N):
            for s in (1,-1): test(c,[(0,1),(j,s)])
        if wmax>=3:
            step=max(1,(N*N)//6000)
            for t,(j,k) in enumerate(itertools.combinations(range(1,N),2)):
                if t%step: continue
                for s1 in (1,-1):
                    for s2 in (1,-1): test(c,[(0,1),(j,s1),(k,s2)])
    logs.sort(); pure.sort(key=lambda t:t[0])
    print("k_{%d^%d}, small split primes %s: %d elements, log10|N|: min %.1f, 1%% quantile %.1f, median %.1f; norms supported on %d and these primes: %d"%(p,n,smalls,cnt,logs[0],logs[len(logs)//100],logs[len(logs)//2],p,len(pure)),flush=True)
    for lg,c,coef,e in pure[:4]: print("    |N| ~ 10^%.1f, c=%d, coef=%s, exponents %s"%(lg,c,coef,e))
run(7,2,[19]); run(43,1,[19]); run(47,1,[53,67,71]); run(71,1,[11]); run(59,1,[53,137])
