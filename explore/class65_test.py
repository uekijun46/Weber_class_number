# Per-prime test for l = 65 mod 128 (class r = 6) in the cyclotomic Z_2-extension.
# If l | h_n/h_{n-1}, then for some prime L | l of Z[zeta_64], every alpha in l*L^{-1} gives a unit eps = (eta_n^alpha)^(1/l) in B_n.
#  height test : |alpha.Log eta_n| / l < sqrt(2^n) * log(2+sqrt 5)                       (bound valid for all relative units)
#  trace test  : alpha(1) odd  and  sum_sigma |eta_n^(sigma alpha)|^(2/l) < 2^n (2^(n+1)-1)  (Komatsu's bound on RE_n^-, Morisawa-Okazaki)
# Either inequality contradicts the existence of eps, so that prime L is excluded.  Floating point: an experiment, not a proof.
import numpy as np, math, sys, time, cypari2
from sympy import nextprime
pari=cypari2.Pari(); pari.allocatemem(10**9,silent=True)
d=32; L0=math.asinh(2)
def v0(n):
    N=2**n; M=2**(n+2); a=np.empty(N,dtype=np.int64); x=1
    for i in range(N): a[i]=x; x=(x*5)%M
    return np.log(np.abs(np.tan(a*math.pi/M)))
LEVELS=range(6,12); V={}; G={}
for n in LEVELS:
    v=v0(n); V[n]=np.array([np.roll(v,-j*2**(n-6)) for j in range(d)]); G[n]=V[n]@V[n].T/2**n
def basis(l,g):
    h=[pow(g,31-i,l) for i in range(d)]
    return np.array([h]+[[l if j==i else 0 for j in range(d)] for i in range(d-1)],dtype=object)
def intlll(B):                                  # exact LLL for the standard inner product (rows of B)
    M=pari("["+";".join(",".join(str(B[j][i]) for j in range(d)) for i in range(d))+"]")   # columns = rows of B
    T=M.qflll()
    T=np.array([[int(T[i][j]) for j in range(d)] for i in range(d)],dtype=object)
    return T@B
CH={}
def lll(B,Gn,n):                                # LLL for the height form: integer LLL on the rounded Cholesky image
    if n not in CH: CH[n]=np.linalg.cholesky(Gn)
    W=np.array(B,dtype=float)@CH[n]*2.0**12
    Wi=[[int(round(x)) for x in row] for row in W]
    M=pari("["+";".join(",".join(str(Wi[j][i]) for j in range(d)) for i in range(d))+"]")
    T=M.qflll()
    T=np.array([[int(T[a][b]) for b in range(d)] for a in range(d)],dtype=object)
    return T@B
def test(l):
    g0=2
    while True:
        g=pow(g0,(l-1)//64,l)
        if pow(g,32,l)==l-1: break
        g0+=1
    res={"height":True,"trace":True}
    for k in range(1,64,2):                  # the 32 primes above l correspond to the 32 primitive 64th roots g^k
        B=intlll(basis(l,pow(g,k,l)))
        for n in LEVELS:
            N=2**n; R=lll(B,G[n],n); Rf=np.array(R,dtype=float)
            ht=min(math.sqrt(r@G[n]@r*N) for r in Rf)/l
            okh= ht<math.sqrt(N)*L0
            cands=list(R)+[R[i]+s*R[j] for i in range(8) for j in range(i+1,12) for s in (1,-1)]
            okt=False
            for a in cands:
                if sum(int(x) for x in a)%2==0: continue
                x=(np.array(a,dtype=float)@V[n])/l
                if np.sum(np.exp(2*x))<N*(2*N-1)*(1-1e-9): okt=True; break
            res["height"]&=okh; res["trace"]&=okt
            if not res["height"] and not res["trace"]: return res
    return res
if __name__=="__main__":
    start=int(float(sys.argv[1])); cnt=int(sys.argv[2]); l=start; done=0; fh=ft=fb=0; t=time.time()
    while done<cnt:
        l=int(nextprime(l))
        if l%128!=65: continue
        r=test(l); done+=1; fh+=not r["height"]; ft+=not r["trace"]; fb+=(not r["height"]) and (not r["trace"])
    print(f"primes l = 65 mod 128 from {start:.3g}: tested {done}; height test fails for {fh}; trace test fails for {ft}; both fail for {fb}   ({time.time()-t:.0f}s)",flush=True)
