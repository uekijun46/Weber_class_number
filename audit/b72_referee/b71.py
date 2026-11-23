import cypari2; pari=cypari2.Pari(); pari.allocatemem(2*10**9)
f=pari('polsubcyclo(49,7)'); print(f)
bnf=pari.bnfinit(f,1); print('h (GRH)',pari('(b)->b.no')(bnf)); print('bnfcertify',pari.bnfcertify(bnf))
print('disc',pari.factor(pari.nfdisc(f)))
# rank minimum
from sympy import primerange, n_order
best=min((l**n_order(l,49),l,n_order(l,49)) for l in primerange(2,2000) if l!=7)
print('min l^f',best)
print([ (l,n_order(l,49)) for l in primerange(2,200) if l!=7 and l**n_order(l,49)<10**5])
