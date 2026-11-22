# Figure-eight knot orbifolds with cone angle 2*pi/2^(n+2):  x = X_n = 2cos(2pi/2^(n+2)),  trace of the longitude  l + 1/l = x^2 - x - 4.
# The longitude eigenvalue l_n is a unit of L_n = B_n(l).  Compute: signature, class number, and whether l_n is a proper power / part of a basis.
import cypari2, math, time
pari=cypari2.Pari(); pari.allocatemem(4*10**9,silent=True)
pari("default(realprecision,60)")
for n in range(0,5):
    N=2**(n+2); t=time.time()
    pari(f"P=x; for(i=1,{n},P=subst(P,x,x^2-2))")           # minimal polynomial of X_n:  X_n^2 - 2 = X_(n-1),  X_0 = 0
    pari("Q=polresultant(subst(P,x,y), l^2-(y^2-y-4)*l+1, y)")  # minimal polynomial of l over Q
    irr=pari("polisirreducible(Q)")
    try:
        pari("alarm(240, b=bnfinit(subst(Q,l,x),1))")
    except Exception as e:
        print(f"n={n}: degree {2**(n+1)}: bnfinit did not finish"); continue
    sig=pari("b.sign"); h=pari("b.no"); rk=pari("#b.fu")
    ev=pari("bnfisunit(b, Mod(x,b.pol))")                    # exponents of l_n on the fundamental units (last entry: torsion)
    g=pari("my(v=bnfisunit(b, Mod(x,b.pol))); gcd(vecextract(v, Str(\"1..\", #v-1)))")
    # lengths: log|l| at the real places
    xs=[2*math.cos(2*math.pi*k/N) for k in range(1,N//2,2)]
    lens=[math.acosh(abs(x*x-x-4)/2) for x in xs if x>-1]
    print(f"n={n}: cone angle 2pi/{N}; deg L_n={2**(n+1)}, irreducible={irr}, signature={sig}, class number={h}, unit rank={rk}; "
          f"gcd of exponents of l_n = {g}; singular length at the orbifold = {lens[0]:.4f}; sum over hyperbolic conjugates = {sum(lens):.4f}  ({time.time()-t:.0f}s)",flush=True)
print("for comparison: 2*vol(figure-eight)*N/(2*pi) with vol=2.029883:  ", [round(2*2.029883*2**(n+2)/(2*math.pi)/2,4) for n in range(5)])
