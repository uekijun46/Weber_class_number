# independent: T = prod_j (x - eta_j), eta_j = sum_{h in g^j H} 2... computed with arb balls; compare to header
from flint import arb, arb_poly, ctx, fmpz_poly
import cypari2; pari=cypari2.Pari()
ctx.prec=800
M=343
H=[x for x in range(1,M) if pow(x,6,M)==1]
g=3  # check primitive root mod 343
assert len({pow(g,k,M) for k in range(294)})==294
reps=[pow(g,j,M) for j in range(49)]
etas=[sum((2*arb.pi()*(r*h%M)/M).cos() for h in H) for r in reps]  # H=-H so period is real: sum of cos
P=arb_poly([1])
for e in etas: P=P*arb_poly([-e,1])
coeffs=[]
for c in P.coeffs():
    m=c.mid(); r=c.rad()
    u=c.unique_fmpz(); n=int(u) if u is not None else None
    assert n is not None, c
    coeffs.append(n)
hdr=[int(x) for x in open('/home/claude/weber_class_number/B72/cert_B72_pilot.txt').readline()[2:].split()]
print('match header:',coeffs==hdr)
T=pari('Pol(Vecrev(%s))'%coeffs)
print('pari polsubcyclo equal:', T==pari('polsubcyclo(343,49)'))
print('irreducible:',pari.polisirreducible(T))
d=pari.poldisc(T); print('v7(disc T)=',pari.valuation(d,7))
r=d/pari('7')**pari.valuation(d,7); print('cofactor is square:',pari.issquare(r), 'log10 cofactor',float(pari.log(pari.abs(r))/pari.log(10)) if r!=1 else 0)
# all headers equal?
import glob
hs={open(f).readline() for f in glob.glob('/home/claude/weber_class_number/B72/cert*.txt')}
print('distinct headers',len(hs))
# field disc via nfdisc with factorization: 
fa=pari.factor(r) if r!=1 else None
print('cofactor factorization', fa)
print('nfdisc', pari.factor(pari.nfdisc(T)))
