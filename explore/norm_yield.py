# How often is the norm of a "small" element of k_m (m = p^n) a product of small split primes?
# Elements: alpha = c + sum a_i eta_i, eta_i = Gaussian periods, a sparse with entries +-1.
import sys, math, random, time
from mpmath import mp, mpf, cos, pi
def setup(p,n):
    M=p**(n+1); N=p**n
    if p==2:
        M=2**(n+2); mp.prec=260; S=200
        E=[int(2*cos(2*pi*pow(5,i,M)/M)*mpf(2)**S) for i in range(N)]
        return M,N,E,S
    g=next(x for x in range(2,M) if x%p and all(pow(x,(p-1)*p**n//r,M)!=1 for r in set([p]+[d for d in range(2,p) if (p-1)%d==0 and all(d%e for e in range(2,d))])))
    H=[pow(g,N*t,M) for t in range(p-1)]                    # subgroup of order p-1
    mp.prec=260
    eta=[sum(cos(2*pi*(pow(g,i,M)*h%M)/M) for h in H) for i in range(N)]
    S=200; E=[int(e*mpf(2)**S) for e in eta]
    return M,N,E,S
def norm(E,S,N,c,idx,sg):
    prod=1
    for j in range(N):
        v=c<<S
        for i,s in zip(idx,sg): v+=s*E[(i+j)%N]
        prod*=v
    sh=S*N; return (prod+(1<<(sh-1)))>>sh if prod>=0 else -((-prod+(1<<(sh-1)))>>sh)
def run(p,n,trials,bounds=(10**4,10**5,10**6),seed=1):
    M,N,E,S=setup(p,n); random.seed(seed)
    lim=max(bounds); sv=bytearray([1])*(lim+1)
    for i in range(2,int(lim**.5)+1):
        if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
    split=[q for q in range(3,lim+1) if sv[q] and q!=p and (pow(q,p-1,M)==1 if p>2 else q%M in (1,M-1))]
    t0=time.time(); logs=[]; smooth={b:0 for b in bounds}; primecof=0; tot=0
    for _ in range(trials):
        w=random.choice((2,3,4)); idx=random.sample(range(N),w); sg=[random.choice((-1,1)) for _ in idx]; c=random.choice((-1,0,1))
        x=abs(norm(E,S,N,c,idx,sg))
        if x==0: continue
        tot+=1; logs.append(math.log10(x))
        while x%p==0: x//=p
        big=0
        for q in split:
            if q*q>x and x<=lim: break
            while x%q==0: x//=q; big=q
            if x==1: break
        if x>1: big=x
        for b in bounds:
            if big<=b: smooth[b]+=1
    logs.sort()
    print("k_{%d^%d}: degree %d, %d elements, %.0f s (%.0f us each); log10|N|: median %.1f, 10%%-90%%: %.1f-%.1f; split primes <= %s: %s"%(p,n,N,tot,time.time()-t0,(time.time()-t0)/tot*1e6,logs[len(logs)//2],logs[len(logs)//10],logs[9*len(logs)//10],bounds,[sum(1 for q in split if q<=b) for b in bounds]))
    print("   fraction with all prime factors <= bound:",{b:"%.2e"%(smooth[b]/tot) for b in bounds},flush=True)
if __name__=="__main__":
    import ast
    for (p,n,t) in ast.literal_eval(sys.argv[1]):
        run(p,n,t)
