#!/usr/bin/env python3
# Audit k_37, part (B): a crude but rigorous bound with mpmath INTERVAL arithmetic (mp.iv), independent of Arb and of relations/.
#   G(c) <= head + sum over subintervals [a,b] of [1e-3, 60] of  sup f([a,b]) (b-a) + tail,
#   where f(x) = (1 - exp(-(x/c)^2))/2 (1/sinh(x/2) + 1/cosh(x/2)) is evaluated on the interval [a,b] (the interval value contains the range),
#   head = (eps^2/2 + eps^3/6)/c^2 for [0, eps], tail = 2 e^(-R/2) (1 + 1/(1 - e^(-R))) for [R, oo)   (derived again in the report).
#   prime sum >= sum_{q in T} sum_{m<=3} 2 log q q^(-m/2) exp(-(m log q / c)^2)   (lower ends of the interval values).
#   B >= pi/2 + gamma + log 8 pi - (72/37) log 37 - G_up + S_low;   h < 2 c sqrt(pi) e^((c/4)^2) / (37 B).
import sys
from mpmath import iv, mp
iv.dps=30
T=[int(l) for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
c=iv.mpf(7); eps=iv.mpf(1)/1000; R=iv.mpf(60)
def f(x):
    u=iv.exp(x/2); v=iv.exp(-x/2)                     # 1/sinh(x/2) = 2/(u - v), 1/cosh(x/2) = 2/(u + v)
    return (1-iv.exp(-(x/c)**2))/2*(2/(u-v)+2/(u+v))
G=(eps**2/2+eps**3/6)/c**2+2*iv.exp(-R/2)*(1+1/(1-iv.exp(-R)))
# subdivision: geometric on [1e-3, 1], uniform on [1, 60]
pts=[iv.mpf(1)/1000*iv.mpf(10)**(iv.mpf(k)/2000) for k in range(0,6001)]
pts=[p_.a for p_ in pts]; pts[0]=mp.mpf(1)/1000; pts[-1]=mp.mpf(1)
pts+= [mp.mpf(1)+mp.mpf(k)*59/20000 for k in range(1,20001)]
for a,b in zip(pts,pts[1:]):
    x=iv.mpf([a,b]); G+=iv.mpf(f(x).b)*(iv.mpf(b)-iv.mpf(a))
S=iv.mpf(0)
for q in T:
    lq=iv.log(q)
    for m in (1,2,3): S+=2*lq/iv.mpf(q)**(iv.mpf(m)/2)*iv.exp(-(m*lq/c)**2)
B=iv.pi/2+iv.euler+iv.log(8*iv.pi)-iv.mpf(72)/37*iv.log(37)-iv.mpf(G.b)+iv.mpf(S.a)
h=2*c*iv.sqrt(iv.pi)*iv.exp((c/4)**2)/(37*iv.mpf(B.a))
print("interval arithmetic (mpmath iv), c = 7:  G(7) <= %s;  prime sum >= %s;  B >= %s;  h < %s"%(mp.nstr(G.b,12),mp.nstr(S.a,12),mp.nstr(B.a,10),mp.nstr(h.b,8)))
