# Untrusted search step for k_p (first layer): generators of the primes above the small split primes, from PARI's class group
# machinery (bnfinit assumes GRH; the output is only a list of candidate rows, to be checked by relations/verify.py).
#   python3 base_pari.py p Y OUT      rows "q | c | i:a ... | q^1" for the split primes q < Y, in the basis 1, eta_0, ..., eta_(p-2)
import sys, time, cypari2
sys.path.insert(0,"../../relations")
from verify import K
p,Y,out=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
F=K(p,1); N=F.N; M=F.M; g=F.g
pari=cypari2.Pari(); pari.allocatemem(4*10**9); pari.set_real_precision(200)
pari("M=%d; N=%d; g=%d; z=exp(2*Pi*I/M); H=vector(%d,t,lift(Mod(g,M)^(N*(t-1))));"%(M,N,g,p-1))
pari("per(i)=real(sum(k=1,#H,z^lift(Mod(g,M)^i*H[k])));")
pari("P=round(real(prod(i=0,N-1,x-per(i))));")
if not bool(pari("polisirreducible(P) && poldegree(P)==N && abs(subst(P,x,per(0)))<1e-100")): raise SystemExit("period polynomial")
pari.set_real_precision(38)
# bnfinit is fast for the reduced polynomial Q of polredbest (75 s), and very slow for P itself
t0=time.time(); pari("QQ=polredbest(P,1); Q=QQ[1]; a0=lift(QQ[2]); bnf=bnfinit(Q);"); print("bnfinit: %.0fs, class group %s"%(time.time()-t0,pari("bnf.cyc")),flush=True)
# Mod(a0,Q) is a root of P, hence a period eta_k; the automorphism s is chosen with s(eta_k) = eta_(k+1)
pari("y=polrootsreal(Q)[1]; v0=subst(a0,x,y); k0=-1; for(k=0,N-1, if(abs(v0-per(k))<1e-15, k0=k)); cj=nfgaloisconj(bnf); s=0; for(k=1,#cj, if(abs(subst(a0,x,subst(cj[k],x,y))-per(k0+1))<1e-15, s=cj[k]));")
if bool(pari("k0<0 || s==0")): raise SystemExit("automorphism not found")
pari("et=vector(N); et[1]=Mod(a0,Q); for(i=2,N, et[i]=subst(lift(et[i-1]),x,Mod(s,Q)));")
pari("A=matrix(N,N,i,j, if(i==1, j==1, polcoeff(lift(et[i-1]),j-1))); Ai=A^(-1);")
pari("gen(q)={my(pr=idealprimedec(bnf,q)[1],r=bnfisprincipal(bnf,pr),a,v,c); if(r[1]!=0*r[1], return([])); a=lift(nfbasistoalg(bnf,r[2])); if(abs(norm(Mod(a,Q)))!=q, return([])); v=vector(N,j,polcoeff(a,j-1)); c=v*Ai; if(denominator(c)!=1, return([])); c};")
f=open(out,"w"); f.write("# field %d 1\n"%p); n=0; t1=time.time()
for q in range(2,Y):
    if not (F.split(q) and bool(pari("isprime(%d)"%q))): continue
    c=pari("gen(%d)"%q)
    if len(c)==0: print("no generator for",q,flush=True); continue
    c=[int(x) for x in c]; coef=" ".join("%d:%d"%(i,c[i+1]) for i in range(N-1) if c[i+1])
    f.write("%d | %d | %s | %d^1\n"%(q,c[0],coef,q)); f.flush(); n+=1
    print("%d: max |coefficient| has %d digits (%.0fs)"%(q,len(str(max(abs(x) for x in c))),time.time()-t1),flush=True)
print("rows: %d"%n)
