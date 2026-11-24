# Floating point exploration (no proof): how far do the constant-factor improvements go?
import numpy as np, math, mpmath as mp, sys
C0=math.pi/2+0.5772156649015329+math.log(8*math.pi)
def gauss(c): return (lambda x: np.exp(-(x/c)**2)), c*math.sqrt(math.pi)/2
def bohman(c):
    def f(x):
        t=np.minimum(np.abs(x)/c,1.0); return (1-t)*np.cos(np.pi*t)+np.sin(np.pi*t)/np.pi
    return f, 4*c/math.pi**2
def G(f):   # integral of (1-F)/2 (1/sinh(x/2)+1/cosh(x/2)), F=f/cosh(x/2)
    x=np.concatenate([np.linspace(1e-9,1,20001),np.linspace(1,60,60001)[1:],np.linspace(60,2000,40001)[1:]])
    F=f(x)/np.cosh(np.minimum(x/2,300)); y=(1-F)/2*(1/np.sinh(np.minimum(x/2,300))+1/np.cosh(np.minimum(x/2,300)))
    return float(np.trapezoid(y,x))
def psum(f,q,N,p,ram):
    L=np.log(q); s=0.0
    for m in (1,2,3): s+=float(np.sum(4*L*f(m*L)/(q**m+1)))     # 2 log q q^(-m/2) F(m log q)
    if ram:
        Lp=math.log(p); s+=sum(4*Lp*float(f(np.array([m*Lp]))[0])/(p**m+1) for m in range(1,40))/N
    return s
def best(q,N,p,logrd,fam,widths,ram):
    out=None
    for c in widths:
        f,Hh=fam(c); B=C0-logrd-Gc[(fam.__name__,c)]+psum(f,q,N,p,ram)
        if B>0:
            h=4*Hh/(N*B)
            if out is None or h<out[0]: out=(h,c,B)
    return out
fields={"k_{3^4}":(3,81,364/81*math.log(3),163,np.load("/home/claude/rel_audit/split_primes.npy").astype(np.float64)),}
# k_{7^2}: split primes <= 2.2e7 by a simple sieve
X7=22000000; sv=np.ones(X7+1,bool); sv[:2]=False
for i in range(2,int(X7**.5)+1):
    if sv[i]: sv[i*i::i]=False
pr=np.nonzero(sv)[0]; res=np.array([pow(int(r),6,343)==1 for r in range(343)]); q7=pr[res[pr%343]].astype(np.float64)
fields["k_{7^2}"]=(7,49,138/49*math.log(7),197,q7)
Wg=[30,40,50,60,70,80,90,100,120,150,200,300,500]; Wb=[60,80,100,120,150,180,220,260,320,400,500,700,1000]
Gc={}
for c in Wg: Gc[("gauss",c)]=G(gauss(c)[0])
for c in Wb: Gc[("bohman",c)]=G(bohman(c)[0])
print("check: G(gauss,100)=%.7f (rigorous value 1.2646927)"%Gc[("gauss",100)])
for name,(p,N,logrd,hmin,q) in fields.items():
    print(name,"split primes:",len(q),"max %.3g"%q[-1])
    for fam,W in ((gauss,Wg),(bohman,Wb)):
        for ram in (False,True):
            b=best(q,N,p,logrd,fam,W,ram)
            # least X (on a grid) such that the primes < X suffice for h < hmin
            lo=None
            for X in np.exp(np.linspace(math.log(q[0])+0.01,math.log(q[-1])+1e-9,60)):
                bb=best(q[q<X],N,p,logrd,fam,W,ram)
                if bb and bb[0]<hmin: lo=(X,bb[1]); break
            # least X such that B > 0 suffices (h < 1e9, if no prime below 1e9 divides h)
            lo2=None
            for X in np.exp(np.linspace(math.log(q[0])+0.01,math.log(q[-1])+1e-9,60)):
                bb=best(q[q<X],N,p,logrd,fam,W,ram)
                if bb and bb[0]<1e9: lo2=(X,bb[1]); break
            print("  %-7s ramified term %-5s: all primes: h < %.2f (width %d);  h < %d already with primes < %.3g (width %s);  B > 0 with primes < %.3g"%(fam.__name__,ram,b[0],b[1],hmin,lo[0] if lo else float('nan'),lo[1] if lo else None,lo2[0] if lo2 else float('nan')),flush=True)
