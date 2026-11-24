# Variant of the bound with one more term (not used in the proofs of the note; see docs/improvements.md).
# The explicit formula behind Miller's lemma contains the sum over ALL prime ideals of the Hilbert class field, and every term
# is positive.  The prime ideal of k_{p^n} above p is principal and has norm p, so it splits completely in the Hilbert class
# field.  Keeping these terms adds
#     R(c) = (2/N) * sum_m log(p) p^(-m/2) F(m log p) = (1/N) * sum_m 4 log(p) exp(-(m log(p)/c)^2) / (p^m + 1)
# to the number B of the lemma.  Ball arithmetic; printed values are rounded in the safe direction.
#   python3 bound_plus.py p n c:S c:S ...       (S = a proven lower bound of the sum over the split primes, as a decimal string)
import sys
from flint import arb, ctx
from verify import K, bound, up, down
def ramified_term(F,c):
    ctx.prec=200; L=arb(F.p).log(); S=arb(0)
    for m in range(1,61): S+=4*L*(-(m*L/arb(c))**2).exp()/(arb(F.p)**m+1)
    return (S/F.N).lower()
if __name__=="__main__":
    p,n=int(sys.argv[1]),int(sys.argv[2]); F=K(p,n); hm=F.hmin(); pi=arb.pi(); logrd,_=F.logrd()
    for item in sys.argv[3:]:
        c,S=item.split(":"); c=int(c); S=arb(S).lower(); I_up,_,_,_=bound(F,[],c); R=ramified_term(F,c)
        out=[]
        for extra in (arb(0),R):
            B=(pi/2+arb.const_euler()+(8*pi).log()-logrd-I_up+S+extra).lower()
            out.append("B >= %s, h < %s"%(down(B,7),up(2*arb(c)*pi.sqrt()/(F.N*B),4)) if B>0 else "B <= 0")
        print("k_{%d^%d}, c=%d: ramified term >= %s; without it: %s; with it: %s   (least possible class number > 1: %d)"%(p,n,c,down(R,7),out[0],out[1],hm))
