# Independent verifier for the certificate of h(B_{7,2}) = 1.  Uses only python-flint (exact integers + ball arithmetic).
# It does NOT trust PARI: the field, every norm, and the analytic bound are recomputed here.
# Usage: python3 verify_cert.py CERT_FILE [CERT_FILE ...]
import sys, time
from flint import fmpz, fmpz_poly, arb, acb, ctx
M=343; DEG=49
Phi=fmpz_poly([1 if k%49==0 else 0 for k in range(295)])   # 1+x^49+...+x^294
H=[x for x in range(1,M) if pow(x,6,M)==1]; assert len(H)==6
E=fmpz_poly([1 if k in H else 0 for k in range(M)]) % Phi   # theta = sum_{h in H} zeta^h  (a Gaussian period)
def mulmod(a,b): return (a*b)%Phi
def eval_at_theta(coeffs):            # Horner in Z[zeta]/(Phi)
    r=fmpz_poly([0])
    for a in reversed(coeffs): r=mulmod(r,E)+fmpz_poly([a])
    return r
def build_T():
    # minimal polynomial of theta: product over the 49 cosets g^j H of (x - theta_j), via power sums would also work;
    # here: characteristic polynomial of multiplication by theta on Z[zeta]/(Phi) is T^6, so take the squarefree part.
    from flint import fmpz_mat
    n=294; rows=[]
    for k in range(n):
        v=mulmod(E,fmpz_poly([0]*k+[1])); c=v.coeffs()+[0]*n; rows.append(c[:n])
    cp=fmpz_mat(rows).charpoly()
    fac=cp.factor()[1]
    assert len(fac)==1 and fac[0][1]==6 and fac[0][0].degree()==DEG, "unexpected factorization"
    return fac[0][0]
def check_rows(fns):
    T=build_T()
    assert eval_at_theta([int(c) for c in T.coeffs()])==0          # T(theta)=0 exactly; T irreducible of degree 49
    primes=[]; t0=time.time()
    import itertools
    for line in itertools.chain.from_iterable(open(fn) for fn in fns):
        if not line.strip(): continue
        if line.startswith("#T"):                                    # optional header: generator's polynomial, low degree first
            assert [int(x) for x in line[2:].split()]==[int(x) for x in T.coeffs()], "generator used a different polynomial"
            continue
        ps,Ds,cs=line.split(";"); p=int(ps); D=int(Ds); a=[int(x) for x in cs.split()]
        assert len(a)==DEG and D>0
        assert fmpz(p).is_prime() and pow(p,6,M)==1, f"{p}: not a split prime"
        A=eval_at_theta(a)                                           # D*alpha as an element of Z[zeta_343]
        co=[int(x) for x in A.coeffs()]
        assert all(x%D==0 for x in co), f"{p}: alpha is not integral"
        A=fmpz_poly([x//D for x in co])                              # alpha in Z[zeta_343], and alpha lies in K
        R=abs(int(Phi.resultant(A)))                                 # = |N_{Q(zeta)/Q}(alpha)| = |N_{K/Q}(alpha)|^6
        assert R==p**6, f"{p}: norm check failed"
        primes.append(p)
    return T,sorted(set(primes)),time.time()-t0
def bound(primes,c):
    ctx.prec=200
    c=arb(c); pi=arb.pi(); gam=arb.const_euler()
    # conductor-discriminant formula: 6 characters of conductor 7^2 and 42 of conductor 7^3,
    # so disc(B_{7,2}) = 7^(6*2+42*3) = 7^138 and the root discriminant is 7^(138/49)
    logrd=arb(138)/49*arb(7).log()
    F=lambda x:(-(x/c)**2).exp()/(x/2).cosh()
    eps=arb(1)/1000; Mx=arb(200)
    f=lambda x,_:(1-(-(x/c)**2).exp()/(x/2).cosh())/2*(1/(x/2).sinh()+1/(x/2).cosh())
    mid=acb.integral(f,acb(eps),acb(Mx)).real
    head=(1/c**2+arb(1)/8)*(eps**2/2+eps**3/6)
    tail=(1/(1-(-Mx).exp())+1)*2*(-Mx/2).exp()
    I_up=(mid+head+tail).upper()
    S=arb(0)
    for p in primes:
        lp=arb(p).log()
        for m in range(1,9): S+=2*lp/arb(p)**(arb(m)/2)*F(m*lp)
    S_low=S.lower()
    B=(pi/2+gam+(8*pi).log()-logrd-I_up+S_low).lower()
    return I_up,S_low,B,( (2*c*pi.sqrt()/(49*B)).upper() if B>0 else None)
if __name__=="__main__":
    T,primes,dt=check_rows(sys.argv[1:])
    print(f"field: degree {T.degree()} subfield of Q(zeta_343), minimal polynomial rebuilt and checked")
    print(f"rows verified: {len(primes)} split primes, largest {primes[-1]}  ({dt:.1f}s)")
    best=None
    for c in (40,50,60,80,100,150):
        I_up,S_low,B,hb=bound(primes,c)
        print(f"c={c}: integral <= {float(I_up):.8f}, prime sum >= {float(S_low):.8f}, B >= {float(B):.6f}"+("" if hb is None else f", h < {float(hb):.2f}"))
        if hb is not None and (best is None or hb<best): best=hb
    if best is None: print("B is not yet positive: more primes are needed")
    elif best<197: print("VERDICT: PASS -- h(B_{7,2}) = 1  (bound below 197, no further input needed)")
    elif best<10000: print("VERDICT: bound below 10^4 -- h = 1 follows once no odd prime < 10^4 other than 7 divides h (Hakkarainen's table, or a direct unit test)")
    else: print("bound not yet below 10^4: more primes are needed")
