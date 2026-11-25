#!/usr/bin/env python3
# Audit k_37, part (B): the numbers of Miller's Lemma 5.2 (under GRH) recomputed with mpmath (50 digits; not interval arithmetic,
# so these are checks of the values, not proofs) and, separately, with python-flint ball arithmetic in a different way
# (the integral by the series  G(c) = int_0^oo (1 - exp(-(x/c)^2)) e^(-x/2) (1/(1-e^(-x)) + 1/(1+e^(-x))) dx  is not used; arb is
# used only for the prime sum with all m until the terms are below 10^-60, giving an enclosure).
# Nothing is imported from relations/.
import sys, mpmath as mp
mp.mp.dps=50
p=37; N=37
T=[int(l) for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
print("primes: %d; smallest %d; largest %d; below 10^6: %d; all distinct: %s; all split (q^36 = 1 mod 1369, q != 37): %s"%(len(T),T[0],T[-1],sum(1 for q in T if q<10**6),len(set(T))==len(T),all(q!=37 and pow(q,36,1369)==1 for q in T)))
logrd=mp.mpf(72)/37*mp.log(37)
const=mp.pi/2+mp.euler+mp.log(8*mp.pi)
print("pi/2 + gamma + log(8 pi) = %s;  log rd = (72/37) log 37 = %s"%(mp.nstr(const,15),mp.nstr(logrd,15)))
def G(c):
    f=lambda x:(1-mp.exp(-(x/c)**2))/2*(1/mp.sinh(x/2)+1/mp.cosh(x/2))
    return mp.quad(f,[0,1,2,4,8,16,32,64,128,256,512])+0   # the tail beyond 512 is below 4 e^(-256)
def S(c,primes,mmax=None):
    s=mp.mpf(0)
    for q in primes:
        lq=mp.log(q); m=1
        while True:
            t=2*lq/mp.mpf(q)**(mp.mpf(m)/2)*mp.exp(-(m*lq/c)**2)
            s+=t; m+=1
            if (mmax and m>mmax) or t<mp.mpf(10)**(-60): break
    return s
def H4(c): return 4*mp.quad(lambda x:mp.exp(-(x/c)**2)*mp.cosh(x/2),[0,c,2*c,4*c,8*c,16*c,mp.inf])
for c in (6,6.5,7,7.5,8):
    c=mp.mpf(c); g=G(c); s=S(c,T); s40=S(c,T,40); B=const-logrd-g+s40
    num=2*c*mp.sqrt(mp.pi)*mp.exp((c/4)**2)
    print("c=%s: G=%s  sum(all m)=%s  sum(m<=40)=%s  B=%s  h < %s   [4 int F cosh(x/2) = %s, 2c sqrt(pi) e^((c/4)^2) = %s]"%(mp.nstr(c,3),mp.nstr(g,13),mp.nstr(s,13),mp.nstr(s40,13),mp.nstr(B,10),mp.nstr(num/(N*B),8),mp.nstr(H4(c),12),mp.nstr(num,12)))
# how much the large primes contribute at c = 7
c=mp.mpf(7)
for X in (10**4,10**5,10**6,10**7,10**10):
    s=S(c,[q for q in T if q<X],40); B=const-logrd-G(c)+s
    print("c=7, primes of T below %d (%d primes): sum = %s, B = %s, h < %s"%(X,sum(1 for q in T if q<X),mp.nstr(s,10),mp.nstr(B,8),mp.nstr(2*c*mp.sqrt(mp.pi)*mp.exp((c/4)**2)/(N*B),6) if B>0 else "-"))
# the split primes of k_37 below 10^6 and how many of them are in T
def isprime(n):
    if n<2: return False
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%a==0: return n==a
    d=n-1; s=0
    while d%2==0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True
H=[r for r in range(1,1369) if r%37 and pow(r,36,1369)==1]
sp=sorted(q for r in H for q in range(r,10**6,1369) if isprime(q)); Ts=set(T)
print("split primes of k_37 below 10^6: %d; in T: %d; the first 20: %s"%(len(sp),sum(1 for q in sp if q in Ts),sp[:20]))
print("smallest split primes NOT in T:",[q for q in sp if q not in Ts][:10])
print("first 16 split primes all in T:",all(q in Ts for q in sp[:16]),"; 16th =",sp[15],"; number below 5100:",sum(1 for q in sp if q<5100))
s16=S(mp.mpf(8),sp[:16],40); B=const-logrd-G(mp.mpf(8))+s16; print("docs/grh.md section 1: c=8, first 16 split primes: h < %s"%mp.nstr(2*8*mp.sqrt(mp.pi)*mp.exp(4)/(N*B),6))
