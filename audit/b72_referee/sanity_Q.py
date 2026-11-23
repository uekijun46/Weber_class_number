# sanity test of the constants on K=Q (n=1, rd=1): explicit formula should give
# pi/2+gamma+log8pi - I(c) + 2 sum_{all p,m} log p p^{-m/2}F(m log p) = 2c sqrt(pi) - sum_rho Phi(rho)
from mpmath import mp, mpf, quad, exp, cosh, sinh, log, sqrt, pi, euler, cos, inf, zetazero
from sympy import primerange
mp.dps=30
for c in [2,3]:
    c=mpf(c)
    F=lambda x: exp(-(x/c)**2)/cosh(x/2)
    I=quad(lambda x:(1-F(x))/2*(1/sinh(x/2)+1/cosh(x/2)),[0,1,5,20,60,inf])
    S=mpf(0)
    for p in primerange(2,3*10**6):
        L=log(p); m=1
        while True:
            t=L*2/(mpf(p)**m+1)*exp(-(m*L/c)**2)
            S+=t; m+=1
            if t<mpf(10)**-25: break
    lhs=pi/2+euler+log(8*pi)-I+2*S
    Z=sum(2*quad(lambda x:2*F(x)*cos(zetazero(k).imag*x),[0,inf]) for k in range(1,4))
    print(float(c),'lhs',lhs,'2c sqrt pi',2*c*sqrt(pi),'diff',lhs-2*c*sqrt(pi),'-zerosum(first 3)',-Z)
