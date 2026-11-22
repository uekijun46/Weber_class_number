import cypari2
pari = cypari2.Pari()
pari.allocatemem(2*10**9)
from sympy import factorint
def Bp(p,n):
    if p==2:
        N=2**(n+2)
        f = pari(f'factor(charpoly(Mod(y+y^({N}-1),polcyclo({N},y)),x))[1,1]')
        return pari.polredabs(f)
    return pari(f'polredabs(polsubcyclo({p**(n+1)},{p**n}))')
def Qm(m):
    pol=None
    for p,v in sorted(factorint(m).items()):
        g=Bp(p,v)
        pol = g if pol is None else pari.polredbest(pari.polcompositum(pol,g)[0])
    if pol is None: return pari('x')
    return pari.polredabs(pol) if pari.poldegree(pol)<=14 else pol
def conductor(m):
    N=1
    for p,v in factorint(m).items():
        N*= 2**(v+2) if p==2 else p**(v+1)
    return N
