# all rows: integrality via coordinates in PARI's maximal order basis (nfinit, certified by nfcertify),
# norm via nfeltnorm. Usage: full_check.py FILE...
import sys, time
import cypari2; pari=cypari2.Pari(); pari.allocatemem(2*10**9)
T=pari('polsubcyclo(343,49)')
nf=pari.nfinit(T)
if len(sys.argv)>1 and sys.argv[1]=='--certify':
    print('nfcertify:',pari.nfcertify(nf), 'disc', pari.factor(pari('(n)->n.disc')(nf))); sys.exit()
t0=time.time(); n=0; bad=[]; ps=[]
for fn in sys.argv[1:]:
    for l in open(fn):
        if not l.strip() or l.startswith('#'): continue
        s,Ds,cs=l.split(';'); p=int(s); D=int(Ds)
        g=pari('Pol(Vecrev([%s]))'%','.join(cs.split()))/D
        v=pari.nfalgtobasis(nf,g)
        integral=int(pari.denominator(v))==1
        N=pari.nfeltnorm(nf,g)
        if not(integral and abs(int(N))==p and pari.isprime(p) and pow(p,6,343)==1): bad.append(p)
        n+=1; ps.append(p)
print(' '.join(sys.argv[1:])[-60:], 'rows',n,'bad',bad,'t',round(time.time()-t0,1), flush=True)
