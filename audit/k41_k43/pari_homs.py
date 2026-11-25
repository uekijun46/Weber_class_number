#!/usr/bin/env python3
# Audit k_41, k_43, another angle: the prime ideals above the primes of T are recomputed with PARI's number field functions, and the
# audit's homomorphisms (rows.py) and patterns J are compared with them.  No homomorphism of rows.py is used to build PARI's objects.
#   nf = nfinit([f, B]):  f = characteristic polynomial of eta_0 (= polsubcyclo(p^2, p)),  B = 1, eta_0, ..., eta_(p-2) written as
#        polynomials g_i in eta_0.  The g_i are computed exactly from the multiplication table, and checked NUMERICALLY (mpmath, 3000 digits,
#        not a proof) against the definition eta_i = sum_{e in E_i} exp(2 pi i e / p^2): |g_i(eta_0) - eta_i| < 10^-1000 for all i.
#        Required: nf.disc = p^(2(p-1)).
#   (1) for every prime q of T: idealprimedec gives p prime ideals of degree 1.  For each of them, the residue map psi: O_K -> F_q is read
#       off from the matrix b of the anti-uniformizer of PARI's prime ideal (x in P iff b x = 0 mod q; rank 1 modulo q; a nonzero row,
#       normalized by psi(1) = 1); for the first primes of each call also from the Hermite normal form (columns q*1 and omega_k - r_k*1).  Required: the p vectors
#       (psi(eta_0), ..., psi(eta_(p-1))) are exactly the p cyclic shifts of the audit's R_q.
#   (2) for the rows that contain q (all of them, or a sample): v_P(alpha) by nfeltval for the p prime ideals P.  Required: the valuation
#       is 1 if the shift of P lies in the audit's J, and 0 otherwise.  (With (1): psi_P = phi tau^j  =>  P = ker(phi tau^j) = s^j(Q_q).)
#   python3 pari_homs.py p A B SECONDS MODE     primes number A..B-1 of the sorted list (or: index B = the primes dividing the index, from number B on;
#        sample B = all primes in a random order with seed B); MODE = homs | val:K (valuations for at most K rows per prime, K = 0: all)
import sys, json, time, random
import cypari2, mpmath as mp
from fieldp import Field, need, read_cert, CERT
pari=cypari2.Pari(); pari.allocatemem(4*10**9,silent=True)
p=int(sys.argv[1]); sel=sys.argv[2]; a0,b0=(0 if sel in ('index','sample') else int(sel)),int(sys.argv[3]); budget=float(sys.argv[4]); mode=sys.argv[5]; F=Field(p); M=F.M; t0=time.time()
# eta_i as polynomials in eta_0
pw=[[1]+[0]*(p-1)]; E0=F.mat[0]
for k in range(1,p): pw.append([sum(pw[-1][s]*E0[s][t] for s in range(p)) for t in range(p)])     # coordinates of eta_0^k
P=pari.matrix(p,p,[x for r in pw for x in r]); need(P.matdet()!=0,"powers of eta_0"); Pinv=P**(-1); den=pari.denominator(Pinv); Pi=Pinv*den
polys=[]
for i in range(p):
    v=F.elem(0,[(i,1)]); polys.append([int(sum(v[t]*Pi[t,k] for t in range(p))) for k in range(p)])    # den * eta_i = sum_k polys[i][k] eta_0^k
mp.mp.dps=3000; deni=int(den)
num=[sum(mp.cos(2*mp.pi*e/M) for e in F.E[i]) for i in range(p)]
need(all(abs(sum(mp.sin(2*mp.pi*e/M) for e in F.E[i]))<mp.mpf(10)**(-2900) for i in range(p)),"eta_i is not real")
x0=num[0]; xp=[mp.mpf(1)]
for k in range(1,p): xp.append(xp[-1]*x0)
err=max(abs(sum(polys[i][k]*xp[k] for k in range(p))/deni-num[i]) for i in range(p)); need(err<mp.mpf(10)**(-1000),"g_i(eta_0) != eta_i numerically")
x=pari("x"); f=pari.matrix(p,p,[c for r in E0 for c in r]).charpoly(); need(f==pari("polsubcyclo(%d,%d)"%(M,p)),"f")
g=[pari.Pol(list(reversed(pl)))/den for pl in polys]; B=[pari(1)]+g[:p-1]
nf=pari.nfinit([f,B]); need(nf[2]==pari(p)**(2*(p-1)),"nfinit: discriminant")
need(pari.nfalgtobasis(nf,pari(1))==pari.Col([1]+[0]*(p-1)),"the first element of nf.zk is not 1")
ecol=[[int(c) for c in pari.nfalgtobasis(nf,g[i])] for i in range(p)]                        # eta_i in PARI's integral basis nf.zk
need(all(sum(ecol[i][k] for i in range(p))==0 for k in range(p)),"sum eta_i != 0 in PARI")
print("k_%d: g_i(eta_0) = eta_i numerically (largest error %s, denominator of %d digits); nfinit: disc = %d^%d; %.0f s"%(p,mp.nstr(err,3),len(str(deni)),p,2*(p-1),time.time()-t0),flush=True)
cert=read_cert(CERT(p),p); rows=json.load(open("out/rows_%d.json"%p)); Rown=json.load(open("out/periods_%d.json"%p))
allq=sorted(int(q) for q in Rown); byq={}
if sel=="index": qs=json.load(open("out/index_primes_%d.json"%p))[b0:]                 # the primes of T that divide the index of Z[eta_0], from number B on
elif sel=="sample": qs=sorted(random.Random(p+b0).sample(allq,min(len(allq),5000)))[:10**9]; random.Random(p*7+b0).shuffle(qs)   # random order, seed B
else: qs=allq[a0:b0]
for n,row in enumerate(rows):
    for q,J in row: byq.setdefault(q,[]).append((n,J))
rng=random.Random(p*1000003+a0); done=0; nval=0; nrow=0
# psi for a prime ideal pr = [q, a, 1, 1, b] of degree 1: b is the matrix of the multiplication by the anti-uniformizer (x in pr iff b x = 0 mod q),
# its rank modulo q is 1 and every nonzero row is a multiple of the residue map; normalized by psi(1) = 1 (the first basis element is 1).
psifun=pari("(pr,q)->my(M=pr[5]*Mod(1,q),r);if(type(pr[5])!=\"t_MAT\"||matrank(M)!=1,return([0]));for(i=1,#M,r=M[i,];if(r!=0,if(r[1]==0,return([0]));return(lift(r/r[1]))));[0]")
for q in qs:
    if time.time()-t0>budget: break
    dec=pari.idealprimedec(nf,q); need(len(dec)==p and all(int(Pr[2])==1 and int(Pr[3])==1 for Pr in dec),"%d: not p prime ideals of degree 1"%q)
    R=Rown[str(q)]; shiftof={tuple(R[j:]+R[:j]):j for j in range(p)}; need(len(shiftof)==p,"%d: shifts"%q); jt=[]
    for t,Pr in enumerate(dec):
        psi=[int(c) for c in psifun(Pr,q)]; need(len(psi)==p and psi[0]==1,"%d: the matrix of the anti-uniformizer does not have rank 1 modulo q"%q)
        if done<3 and t<3:                                                                    # the same map from the Hermite normal form
            Hn=pari.idealhnf(nf,Pr); need(int(Hn[0,0])==q and all(int(Hn[k,k])==1 for k in range(1,p)),"%d: Hermite normal form"%q)
            need(psi[1:]==[(-int(Hn[0,k]))%q for k in range(1,p)],"%d: the two ways to read off the residue map differ"%q)
        need(sum(c*s_ for c,s_ in zip([int(c) for c in Pr[1]],psi))%q==0,"%d: psi does not vanish on the second generator"%q)
        vec=tuple(sum(c*s_ for c,s_ in zip(ecol[i],psi))%q for i in range(p))
        need(vec in shiftof,"%d: PARI's residue map is not a cyclic shift of the audit's homomorphism"%q); jt.append(shiftof[vec])
    need(sorted(jt)==list(range(p)),"%d: the p prime ideals do not give the p shifts"%q)
    if mode.startswith("val"):
        K=int(mode.split(":")[1]); lst=byq[q] if K==0 or len(byq[q])<=K else rng.sample(byq[q],K)
        for n,J in lst:
            c,coef,fac=cert[n]; col=[0]*p; col[0]=c
            for i,a in coef:
                for k in range(p): col[k]+=a*ecol[i][k]
            al=pari.Col(col); V=[int(pari.nfeltval(nf,al,Pr)) for Pr in dec]
            need(all(V[t]==(1 if jt[t] in J else 0) for t in range(p)) and sum(V)==dict(fac)[q],"%d, line %d: valuations %s"%(q,n,V)); nval+=p; nrow+=1
    done+=1
print("k_%d, selection %s %s (number %d..%d in it; first q = %d, last q = %d), mode %s: %d primes: PARI's %d residue maps are exactly the cyclic shifts of the audit's homomorphism%s; %.0f s"%(p,sel,sys.argv[3],a0,a0+done-1,qs[0],qs[done-1] if done else 0,mode,done,p,"" if not mode.startswith("val") else "; %d (row, prime) entries, %d valuations by nfeltval, all 0/1 and equal to the audit's patterns"%(nrow,nval),time.time()-t0),flush=True)
