import sys,time,json
from fields import *
from miller import *
mp.mp.dps=20
m=int(sys.argv[1]); Qmax=int(sys.argv[2])
P=Qm(m); bnf=pari.bnfinit(pari.nfinit([P,[p for p in factorint(m)]]),1); nf=bnf
split=[]
t=time.time()
for q in primerange(2,Qmax+1):
    if q in factorint(m): continue
    if frob_order_lcm(q,m)!=1: continue
    pr=pari.idealprimedec(nf,q)[0]
    r=pari.bnfisprincipal(bnf,pr,1)
    assert pari.length(r[0]) == 0 or all(int(x)==0 for x in r[0])
    a=r[1]
    # exact checks: integral and norm +-q
    assert all(pari.denominator(x)==1 for x in a), 'not integral'
    N=pari.nfeltnorm(nf,a)
    assert abs(int(N))==q, (q,N)
    split.append(q)
json.dump(split,open(f'split{m}.json','w'))
print(m,'certified',len(split),'split primes <=',Qmax,'time',time.time()-t)
