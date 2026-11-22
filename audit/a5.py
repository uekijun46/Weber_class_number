from fields import *
pari('''fermatmat(bnf,p)=
 my(g=bnf.pol, fu=bnf.fu, U=apply(lift,fu), vd, prec, F, gv, M);
 vd=vecmax(concat([0],apply(u->valuation(denominator(content(u)),p),U)));
 prec=6+poldegree(g)*vd;
 F=factorpadic(g,p,prec);
 if(vecmax(F[,2])>1, error("ramified"));
 gv=vector(#F[,1],i,Pol(apply(c->truncate(c),Vec(F[i,1]))));
 M=matrix(#fu,#gv,i,j, my(u=U[i], N);
    N=polresultant(gv[j],u); N=Mod(numerator(N),p^(prec))/Mod(denominator(N),p^prec) ;
    N=lift(N)%(p^2); if(N%p==0,error("nonunit"));
    lift(Mod((lift(Mod(N,p^3)^(p-1))-1)/p,p)));
 [#gv, if(#fu,matrank(M*Mod(1,p)),0), M, apply(poldegree,gv)]''')
def bnfK(m):
    P=Qm(m)
    return pari.bnfinit(pari.nfinit([P,[p for p in factorint(m)]]),1)
def fm(bnf,p): return pari('fermatmat')(bnf,p)
