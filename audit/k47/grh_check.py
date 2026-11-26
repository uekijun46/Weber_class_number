#!/usr/bin/env python3
# Audit k_47 (copied from audit/k41_k43/grh_check.py; c = 7 replaced by c = 6.5; one line added for T without two primes).
# Audit k_41, k_43, part (B): the numbers of Miller's Lemma 5.2 (under GRH) recomputed without relations/ and without Arb.
#  (1) mpmath, 50 digits (NOT interval arithmetic: a check of the values, not a proof):
#        G(c) = int_0^oo (1 - e^(-(x/c)^2))/2 (1/sinh(x/2) + 1/cosh(x/2)) dx   in two ways: tanh-sinh quadrature on subintervals, and
#        Gauss-Legendre quadrature after the substitution x = e^u is NOT used; instead the second way is the series
#        (1/2)(1/sinh(x/2) + 1/cosh(x/2)) = sum_{k>=0} (e^(-(2k+1)x/2) + (-1)^k e^(-(2k+1)x/2)) = 2 sum_{k even} e^(-(2k+1)x/2),
#        int_0^oo (1 - e^(-(x/c)^2)) e^(-a x) dx = 1/a - (c sqrt(pi)/2) e^(a^2 c^2/4) erfc(a c/2),   a = (4j+1)/2, j >= 0.
#        prime sum with all m (until the terms are below 10^-60) and with m <= 40.
#  (2) mpmath interval arithmetic (mp.iv), crude but rigorous enclosures: upper Riemann sums for G on [1e-3, 60] with the head and
#        tail bounds derived in docs/audit_k37.md section 1.7, the prime sum with m <= 12 (a lower bound, all terms are positive).
#   python3 grh_check.py p
import sys, mpmath as mp
from mpmath import iv
from fieldp import need, is_prime_mr, PRIMES
p=int(sys.argv[1]); N=p; M=p*p; mp.mp.dps=50
lines=[l.strip() for l in open(PRIMES(p))]
T=[int(l) for l in lines if l and not l.startswith("#")]
print("k_%d: primes in span_primes.txt: %d; smallest %s; largest %d; below 10^6: %d; sorted and distinct: %s; all prime (Miller-Rabin, deterministic): %s; all split (q != p, q^(p-1) = 1 mod p^2): %s"%(p,len(T),T[:6],T[-1],sum(1 for q in T if q<10**6),T==sorted(set(T)),all(is_prime_mr(q) for q in T),all(q!=p and pow(q,p-1,M)==1 for q in T)))
need(T==sorted(set(T)) and all(is_prime_mr(q) and q!=p and pow(q,p-1,M)==1 for q in T),"the list of primes")
logrd=mp.mpf(2*(p-1))/p*mp.log(p); const=mp.pi/2+mp.euler+mp.log(8*mp.pi)
print("pi/2 + gamma + log(8 pi) = %s;  log rd = (%d/%d) log %d = %s"%(mp.nstr(const,15),2*(p-1),p,p,mp.nstr(logrd,15)))
def G(c):
    f=lambda x:(1-mp.exp(-(x/c)**2))/2*(1/mp.sinh(x/2)+1/mp.cosh(x/2))
    return mp.quad(f,[0,1,2,4,8,16,32,64,128,256,512])
def G_series(c):
    s=mp.mpf(0); j=0
    while True:
        a=mp.mpf(4*j+1)/2; t=1/a-c*mp.sqrt(mp.pi)/2*mp.exp(a*a*c*c/4)*mp.erfc(a*c/2)      # ~ 2/(a^3 c^2) for large a
        s+=2*t; j+=1
        if j>200000: break
    # tail of sum_{j>J} 2 * 2/(a^3 c^2) is below 4/c^2 * int_J^oo (2j)^(-3) dj; reported with the value
    return s, 4/c**2/(16*mp.mpf(j)**2)
def S(c,primes,mmax=None):
    s=mp.mpf(0)
    for q in primes:
        lq=mp.log(q); m=1
        while True:
            t=2*lq/mp.mpf(q)**(mp.mpf(m)/2)*mp.exp(-(m*lq/c)**2); s+=t; m+=1
            if (mmax and m>mmax) or t<mp.mpf(10)**(-60): break
    return s
best=None
for c2 in range(8,29):
    c=mp.mpf(c2)/2; g=G(c); s40=S(c,T,40); B=const-logrd-g+s40; num=2*c*mp.sqrt(mp.pi)*mp.exp((c/4)**2); h=num/(N*B) if B>0 else None
    print("  c=%4s: G=%s  sum(m<=40)=%s  B=%s  h < %s"%(mp.nstr(c,3),mp.nstr(g,13),mp.nstr(s40,13),mp.nstr(B,10),mp.nstr(h,8) if h else "-"))
    if h and (best is None or h<best[0]): best=(h,c)
print("mpmath: the smallest bound is at c = %s: h < %s"%(mp.nstr(best[1],3),mp.nstr(best[0],10)))
c=mp.mpf(13)/2; g=G(c); gs,gt=G_series(c); s=S(c,T); s40=S(c,T,40); B=const-logrd-g+s40
H4=4*mp.quad(lambda x:mp.exp(-(x/c)**2)*mp.cosh(x/2),[0,c,2*c,4*c,8*c,16*c,mp.inf]); num=2*c*mp.sqrt(mp.pi)*mp.exp((c/4)**2)
print("c=6.5: G (quadrature) = %s;  G (series, 200001 terms) = %s + (tail below %s);  sum(all m) = %s;  sum(m<=40) = %s"%(mp.nstr(g,15),mp.nstr(gs,15),mp.nstr(gt,3),mp.nstr(s,15),mp.nstr(s40,15)))
print("c=6.5: B = %s;  2c sqrt(pi) e^((c/4)^2) = %s  (4 int_0^oo F(x) cosh(x/2) dx = %s);  h < %s"%(mp.nstr(B,12),mp.nstr(num,12),mp.nstr(H4,12),mp.nstr(num/(N*B),10)))
for X in (10**4,10**5,10**6,10**7,10**10):
    sub=[q for q in T if q<X]; sX=S(c,sub,40); BX=const-logrd-g+sX
    print("c=6.5, primes of T below %d (%d primes): sum = %s, B = %s, h < %s"%(X,len(sub),mp.nstr(sX,10),mp.nstr(BX,8),mp.nstr(num/(N*BX),6) if BX>0 else "-"))
Tp=[q for q in T if q not in (6453089,184398551)]; sp_=S(c,Tp,40); Bp=const-logrd-g+sp_
print("c=6.5, T without 6453089 and 184398551 (%d primes): sum = %s, B = %s, h < %s"%(len(Tp),mp.nstr(sp_,15),mp.nstr(Bp,12),mp.nstr(num/(N*Bp),10)))
# (2) interval arithmetic
iv.dps=30; ci=iv.mpf(13)/2; eps=iv.mpf(1)/1000; Rr=iv.mpf(60)
def f(x):
    u=iv.exp(x/2); v=iv.exp(-x/2); return (1-iv.exp(-(x/ci)**2))/2*(2/(u-v)+2/(u+v))
Gi=(eps**2/2+eps**3/6)/ci**2+2*iv.exp(-Rr/2)*(1+1/(1-iv.exp(-Rr)))
pts=[(iv.mpf(1)/1000*iv.mpf(10)**(iv.mpf(k)/2000)).a for k in range(0,6001)]; pts[0]=mp.mpf(1)/1000; pts[-1]=mp.mpf(1)
pts+=[mp.mpf(1)+mp.mpf(k)*59/20000 for k in range(1,20001)]
need(all(a<b for a,b in zip(pts,pts[1:])) and pts[0]==mp.mpf(1)/1000 and pts[-1]==60,"subdivision")
for a,b in zip(pts,pts[1:]):
    x=iv.mpf([a,b]); Gi+=iv.mpf(f(x).b)*(iv.mpf(b)-iv.mpf(a))
Si=iv.mpf(0)
for q in T:
    lq=iv.log(q)
    for m in range(1,13): Si+=2*lq/iv.mpf(q)**(iv.mpf(m)/2)*iv.exp(-(m*lq/ci)**2)
Bi=iv.pi/2+iv.euler+iv.log(8*iv.pi)-iv.mpf(2*(p-1))/p*iv.log(p)-iv.mpf(Gi.b)+iv.mpf(Si.a)
hi=2*ci*iv.sqrt(iv.pi)*iv.exp((ci/4)**2)/(N*iv.mpf(Bi.a))
print("interval arithmetic (mpmath iv), c = 6.5:  G(6.5) <= %s;  prime sum >= %s;  B >= %s;  h < %s"%(mp.nstr(mp.mpf(Gi.b),12),mp.nstr(mp.mpf(Si.a),12),mp.nstr(mp.mpf(Bi.a),10),mp.nstr(mp.mpf(hi.b),8)))
# the split primes below 10^6 and how many of them are in T; the claims of docs/grh.md section 1
Hs=[r for r in range(1,M) if r%p and pow(r,p-1,M)==1]
sp=sorted(q for r in Hs for q in range(r,10**6,M) if is_prime_mr(q)); Ts=set(T)
print("split primes of k_%d below 10^6: %d; in T: %d; the first 16: %s"%(p,len(sp),sum(1 for q in sp if q in Ts),sp[:16]))
print("smallest split primes NOT in T:",[q for q in sp if q not in Ts][:8])
mp.mp.dps=30
for cnt,cc in ((13,None),(8,None),(16,None)):
    bb=None
    for c2 in range(8,40):
        c_=mp.mpf(c2)/2; B_=const-logrd-G(c_)+S(c_,sp[:cnt],40)
        if B_>0:
            h_=2*c_*mp.sqrt(mp.pi)*mp.exp((c_/4)**2)/(N*B_)
            if bb is None or h_<bb[0]: bb=(h_,c_)
    print("docs/grh.md section 1: the first %d split primes (%d..%d): best h < %s at c = %s"%(cnt,sp[0],sp[cnt-1],mp.nstr(bb[0],6) if bb else "-",mp.nstr(bb[1],3) if bb else "-"))
