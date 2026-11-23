# Rigorous enclosure of I(c)=int_0^inf (1-F)/2 (1/sinh(x/2)+1/cosh(x/2)) dx, F=exp(-(x/c)^2)/cosh(x/2)
# Method: Taylor models with Lagrange remainder via arb_series (not acb.integral).
import sys
from flint import arb, arb_series, ctx
ctx.prec=128
K=14
def series_f(x0,c):
    ctx.cap=K+1
    t=arb_series([x0,1])
    e=(t/2).exp(); ei=1/e
    sh=(e-ei)/2; ch=(e+ei)/2
    F=(-(t/c)*(t/c)).exp()/ch
    return ((1-F)/2*(1/sh+1/ch)).coeffs()
def enclose(c, a=arb(1)/10000, Mx=60, h=arb(1)/4):
    c=arb(c)
    # head [0,a]: 0 <= f(x) <= (1/c^2+1/8)(x+x^2/2)  (proved by hand, see report)
    head_hi=(1/c**2+arb(1)/8)*(a**2/2+a**3/6)
    # tail [Mx,inf): 0 <= f <= (1/2)(1/sinh(x/2)+1/cosh(x/2)) <= (1/(1-e^-Mx)+1) e^{-x/2}
    tail_hi=(1/(1-(-arb(Mx)).exp())+1)*2*(-arb(Mx)/2).exp()
    tot=arb(0); x=a
    while x<Mx:
        b=min(x+min(h,x/2),arb(Mx))
        x0=(x+b)/2; r=(b-x)/2
        cs=series_f(x0,c)
        ball=x.union(b)
        cK=series_f(ball,c)[K]
        s=sum(2*cs[j]*r**(j+1)/(j+1) for j in range(0,K,2))
        rem=cK*2*r**(K+1)/(K+1)
        tot+= s + arb(0,abs(rem).upper())
        x=b
    lo=tot.lower(); hi=(tot+head_hi+tail_hi).upper()
    return lo,hi
if __name__=='__main__':
    for c in map(int,sys.argv[1:]):
        lo,hi=enclose(c); print(c, lo.str(20,radius=False), hi.str(20,radius=False))
