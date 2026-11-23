# independent re-check of certificate rows: alpha=(1/D)sum a_i theta^i in Q[x]/(T)
# integrality: characteristic polynomial of alpha over Q has integer coefficients (=> alpha in O_K)
# norm: N_{K/Q}(alpha) = (-1)^49 * charpoly(0), require |.|=p.  Also p prime (PARI isprime = proven), p^6=1 mod 343.
import sys, random, glob, time
import cypari2; pari=cypari2.Pari(); pari.allocatemem(10**9)
files=sorted(glob.glob('/home/claude/weber_class_number/B72/cert*.txt'))
idx=[]  # (file, line number)
for fn in files:
    for i,l in enumerate(open(fn)):
        if l.strip() and not l.startswith('#'): idx.append((fn,i))
seed=int(sys.argv[1]); n=int(sys.argv[2]); start=int(sys.argv[3]); stop=int(sys.argv[4])
random.seed(seed); samp=random.sample(idx,n)[start:stop]
want={}
for fn,i in samp: want.setdefault(fn,set()).add(i)
T=pari('polsubcyclo(343,49)')  # checked separately against an arb product of Gaussian periods
ok=0; t0=time.time(); bad=[]
for fn,S in want.items():
    for i,l in enumerate(open(fn)):
        if i not in S: continue
        ps,Ds,cs=l.split(';'); p=int(ps); D=int(Ds); a=[int(x) for x in cs.split()]
        assert len(a)==49
        g=pari('Pol(Vecrev(%s))'%('['+','.join(cs.split())+']'))
        al=pari.Mod(g/D,T)
        cp=pari.charpoly(al)
        integral = all(cf.type()=='t_INT' for cf in pari.Vec(cp)) and int(pari.denominator(cp))==1
        N=-pari.polcoef(cp,0)
        good = integral and abs(int(N))==p and bool(pari.isprime(p)) and pow(p,6,343)==1 and pari.poldegree(cp)==49
        # extra: minimal polynomial degree (alpha generates K or a subfield; irrelevant but informative)
        if good: ok+=1
        else: bad.append((fn,i,p,integral,int(N)))
print('checked',len(samp),'ok',ok,'bad',bad,'time',round(time.time()-t0,1))
print('primes', sorted(int(l.split(';')[0]) for fn,S in want.items() for i,l in enumerate(open(fn)) if i in S)[:5], '...')
