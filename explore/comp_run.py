# Class number 1 for composite layers K = prod B_{p_i,n_i} of the cyclotomic \hat{Z}-extension, by induction on subfields.
import cypari2, time, sys, math, numpy as np
import compos2
from compos2 import setup, logrd, B, order, PR
pari=cypari2.Pari(); pari.allocatemem(3*10**9); pari.default("realprecision",150)
def layer(p,n):
    if p==2:
        P=pari("x")
        for _ in range(n): P=pari.subst(P,"x",pari("x^2-2"))
        return P
    P=pari('polsubcyclo(%d,%d)'%(p**(n+1),p**n))
    return P[0] if pari.type(P)=='t_VEC' else P
def field(F):
    T=layer(*F[0])
    for p,n in F[1:]:
        T=pari.polcompositum(T,layer(p,n))
        assert len(T)==1; T=T[0]
    return pari.polredbest(T)
def genus_ok(F,p):          # p-part for a ramified p (layer n=1) that splits completely in the other factors: Chevalley + Hasse
    Fo=[(q,m) for q,m in F if q!=p]; T=field(Fo); t=int(pari.poldegree(T))
    bnf=pari.bnfinit(T,1); fac=pari.factorpadic(T,p,8)[0]; nv=len(fac)
    M=[]
    for u in list(bnf.bnf_get_fu())+[bnf.bnf_get_tu()[1]]:
        row=[]
        for gv in fac:
            v=int(pari.lift(pari.polresultant(gv,pari.lift(u))+pari("O(%d^3)"%p)))%p**2     # local norm of u at the place gv
            assert v%p
            row.append(((pow(v,p-1,p*p)-1)//p)%p)
        M.append(row)
    t=nv
    rk=int(pari.matrank(pari.matrix(len(M),t,[x for row in M for x in row])*pari("Mod(1,%d)"%p)))
    return rk==t-1, rk, t
cs=(6,8,10,12,15,20,25,30,40,50,60,80,100)
def hb(F,X,S):
    N=S[0]; best=math.inf
    for c in cs:
        b=B(F,c,X,S)
        if b>0: best=min(best,2*c*math.sqrt(math.pi)/(N*b))
    return best
def run(F,limit):
    t0=time.time(); S=setup(F); N,qs,g,f,rows,hm,wief=S
    hmin=min(int(l)**order(int(l),N) for l in PR[:400] if N%int(l))
    note=""
    for p,gg in wief:
        ok,rk,t=genus_ok(F,p); note+=" [%d-part: %d primes above %d, unit rank mod %d = %d -> %s]"%(p,t,p,p,rk,"trivial" if ok else "NONTRIVIAL")
        if not ok: return "%s N=%d: %s"%(F,N,note)
    split=qs[g==N]; X=0; ncert=0; h=hb(F,0,S)
    if h>=2 and h>=hmin or (h>=hmin):
        T=field(F); bnf=pari.bnfinit(T,1); i=0
        while h>=hmin:
            if time.time()-t0>limit or i>=len(split): return "%s N=%d rd=%.0f: NOT finished, X=%d certs=%d bound=%.1f hmin=%d%s"%(F,N,math.exp(logrd(F)),X,ncert,h,hmin,note)
            for _ in range(max(5,ncert//3)):
                if i>=len(split): break
                q=int(split[i]); i+=1
                P=pari.idealprimedec(bnf,q); assert len(P)==N
                r=pari.bnfisprincipal(bnf,P[0])
                if any(int(e) for e in r[0]): return "%s: NON-PRINCIPAL prime above %d"%(F,q)
                assert abs(int(pari.norm(pari.nfbasistoalg(bnf,r[1]))))==q
                ncert+=1; X=q
            h=hb(F,X,S)
    return "%s N=%d rd=%.0f: h < %.2f < hmin=%d with %d certificates (split primes <= %d)%s  [%.0fs]"%(F,N,math.exp(logrd(F)),h,hmin,ncert,X,note,time.time()-t0)
if __name__=="__main__":
    Fs=[[(2,1),(3,1)],[(2,2),(3,1)],[(2,1),(5,1)],[(3,1),(5,1)],[(2,1),(7,1)],[(2,1),(3,2)],[(2,2),(5,1)],[(3,1),(7,1)],[(2,3),(3,1)],
        [(2,1),(11,1)],[(2,2),(7,1)],[(2,2),(3,2)]]
    if len(sys.argv)>1: Fs=eval(sys.argv[1])
    for F in Fs: print(run(F,float(sys.argv[2]) if len(sys.argv)>2 else 60),flush=True)
