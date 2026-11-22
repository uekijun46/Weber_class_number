import numpy as np, math
from scipy import integrate, optimize
import feas
from feas import P, split_primes, rd, C0
L=feas.L
def Ith(c,th):
    def f(x):
        if x>700: return 0.0
        F=math.exp(-(x/c)**2)/math.cosh(th*x) if th>0 else math.exp(-(x/c)**2)
        return (1-F)/2*(1/math.sinh(x/2)+1/math.cosh(x/2))
    return integrate.quad(f,1e-12,60,limit=500)[0]+2*math.exp(-30)
def pole(c,th):
    if th>=0.5: return 2*c*math.sqrt(math.pi)
    def f(x):
        a=-(x/c)**2+x/2-th*x
        return 4*math.exp(a)*(1+math.exp(-x))/(1+math.exp(-2*th*x)) if a>-700 else 0.0
    return integrate.quad(f,0,40*c+50,limit=800)[0]
def wsum(sp,c,th,X,N,p0):
    q=sp[sp<=min(X,L)].astype(float); t=np.log(q); s=0.0
    for m in range(1,8):
        s+=np.sum(2*t/q**(m/2)*np.exp(-(m*t/c)**2)/np.cosh(th*m*t))
    if X>L:
        s+=integrate.quad(lambda t:2*math.exp(t/2-(t/c)**2)/math.cosh(th*t)/N,math.log(L),math.log(X),limit=300)[0]
    for m in range(1,60):  # ramified prime
        t=m*math.log(p0); s+=2/N*math.log(p0)/p0**(m/2)*math.exp(-(t/c)**2)/math.cosh(th*t)
    return s
def hbound(p0,n,X,th,sp,cs):
    N=p0**n; r=rd(p0,n); best=math.inf
    for c in cs:
        B=C0-math.log(r)-Ith(c,th)+wsum(sp,c,th,X,N,p0)
        if B>0: best=min(best,pole(c,th)/(N*B))
    return best
def minX(p0,n,th,target,sp,cs,lo=3.0,hi=40.0):
    f=lambda lx: hbound(p0,n,10**lx,th,sp,cs)
    if f(hi)>target: return None
    for _ in range(22):
        mid=(lo+hi)/2
        if f(mid)<target: hi=mid
        else: lo=mid
    return hi
cs_unc=(20,30,40,50,60,80,100,150,200,300,500,1000)
print("== (1) unconditional: log10 of X needed so that h < target ==")
for (p0,n,hm) in [(7,2,197),(29,1,59),(3,4,163),(31,1,311)]:
    sp=split_primes(p0,n); N=p0**n
    row=[]
    for tgt in (hm,1e4,1e6,1e9):
        x=minX(p0,n,0.5,tgt,sp,cs_unc); 
        cnt=(np.sum(sp<=10**x) if 10**x<=L else 10**x/(x*math.log(10))/N)
        row.append("h<%g: X=1e%.2f (#%.2e)"%(tgt,x,cnt))
    print(f"B_({p0},{n}):  "+" | ".join(row))
print("== (2) if zeros of zeta_H satisfy |Re-1/2|<=theta: log10 X needed for h<hmin ==")
for (p0,n,hm) in [(7,2,197),(29,1,59),(3,4,163),(2,7,1e9),(5,3,251),(2,8,1e9)]:
    sp=split_primes(p0,n); N=p0**n; out=[]
    for th in (0.5,0.45,0.4,0.3,0.2,0.1,0.0):
        cs=(20,30,50,80,150,300) if th>=0.45 else np.arange(3,40,1.0)
        x=minX(p0,n,th,hm,sp,cs,lo=1.0,hi=45.0)
        out.append("th=%.2f:%s"%(th,"-" if x is None else "1e%.1f"%x))
    print(f"B_({p0},{n}) first split {sp[:3].tolist()}:  "+"  ".join(out))
