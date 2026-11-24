# Upper bound of the class number UNDER THE GENERALIZED RIEMANN HYPOTHESIS (Miller, Acta Arith. 164 (2014), Lemma 5.2).
# Results of this program are conditional; they are kept apart from the unconditional results.
#   F(x) = exp(-(x/c)^2),  G(c) = int_0^oo (1-F(x))/2 (1/sinh(x/2) + 1/cosh(x/2)) dx,
#   B = pi/2 + gamma + log(8 pi) - log rd - G(c) + 2 sum_{q in S} sum_m log(q) q^(-m/2) F(m log q),
#   B > 0  =>  h < 2 c sqrt(pi) exp((c/4)^2) / (N B)       (S: primes that split completely into principal prime ideals)
# Ball arithmetic.  G(c) is bounded from above as in the unconditional case: on [0, eps] by (eps^2/2 + eps^3/6)/c^2
# (from 1 - F(x) <= x^2/c^2, 1/sinh(x/2) <= 2/x, 1/cosh(x/2) <= 1), on [R, oo) by 2 exp(-R/2) (1 + 1/(1 - exp(-R))).
# The sum is a partial sum (m <= 40) of positive terms.
#   python3 bound_grh.py p n FILE [number of primes]     FILE: rows "q | ..." that were verified by verify.py --complete
import sys
from flint import arb, acb, ctx
from verify import K, up, down, need
def bound_grh(F,primes,c):
    ctx.prec=200; c=arb(c); pi=arb.pi(); logrd,_=F.logrd(); eps=arb(1)/1000; R=arb(200)
    f=lambda x,_:(1-(-(x/c)**2).exp())/2*(1/(x/2).sinh()+1/(x/2).cosh())
    G_up=(acb.integral(f,acb(eps),acb(R)).real+(eps**2/2+eps**3/6)/c**2+(1/(1-(-R).exp())+1)*2*(-R/2).exp()).upper()
    S=arb(0)
    for q in primes:
        lq=arb(q).log()
        for m in range(1,41): S+=2*lq/arb(q)**(arb(m)/2)*(-(m*lq/c)**2).exp()
    B=(pi/2+arb.const_euler()+(8*pi).log()-logrd-G_up+S.lower()).lower()
    return G_up,S.lower(),B,((2*c*pi.sqrt()*((c/4)**2).exp()/(F.N*B)).upper() if B>0 else None)
if __name__=="__main__":
    p,n,fn=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]; F=K(p,n)
    qs=sorted(int(l.split("|")[0]) for l in open(fn) if l.strip() and not l.startswith("#") and not l.startswith("L"))
    need(all(F.split(q) for q in qs) and len(set(qs))==len(qs),"the file must list distinct split primes")
    if len(sys.argv)>4: qs=qs[:int(sys.argv[4])]
    hm=F.hmin(); best=None
    print("UNDER GRH.  field k_{%d^%d}; %d split primes from %s (smallest %d, largest %d); least possible class number > 1: %d"%(p,n,len(qs),fn,qs[0],qs[-1],hm))
    for c2 in range(8,29):                      # c = 4, 4.5, ..., 14
        c=arb(c2)/2; G_up,S_low,B,hb=bound_grh(F,qs,c)
        if hb is not None and (best is None or hb<best[0]): best=(hb,c2/2,G_up,S_low,B)
    need(best is not None,"B <= 0 for every c")
    hb,c,G_up,S_low,B=best
    print("c=%s: G(c) <= %s, prime sum >= %s, B >= %s, h < %s"%(c,up(G_up,10),down(S_low,10),down(B,7),up(hb,4)))
    print("VERDICT (under GRH): PASS  (h < %s < %d)"%(up(hb,4),hm) if hb<hm else "VERDICT (under GRH): bound not below %d"%hm)
