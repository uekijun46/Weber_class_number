import cypari2, time, sys
p=int(sys.argv[1]); K=int(sys.argv[2])
pari=cypari2.Pari(); pari.allocatemem(4*10**9,silent=True); pari.default("realprecision",150)
t0=time.time(); T=pari("polsubcyclo(%d,%d)"%(p*p,p)); 
print("k_%d: disc = %d^%d (10^%.1f)"%(p,p,2*(p-1),2*(p-1)*__import__('math').log10(p)),flush=True)
bnf=pari.bnfinit(T,1); print("bnfinit %.0fs, class group (GRH):"%(time.time()-t0), bnf.bnf_get_cyc(),flush=True)
out=open("cert_k%d.txt"%p,"w"); out.write("#T "+" ".join(str(c) for c in pari.Vecrev(T))+"\n")
q=1; n=0
while n<K:
    q=int(pari.nextprime(q+1))
    if q==p or pow(q,p-1,p*p)!=1: continue
    t=time.time(); P=pari.idealprimedec(bnf,q); r=pari.bnfisprincipal(bnf,P[0])
    if any(int(e) for e in r[0]): print("NON-PRINCIPAL?",q); n+=1; continue
    a=pari.lift(pari.nfbasistoalg(bnf,r[1])); D=int(pari.denominator(pari.content(a))); A=a*D
    assert abs(int(pari.polresultant(T,A)))==q*D**p
    out.write("%d ; %d ; %s\n"%(q,D," ".join(str(c) for c in pari.Vecrev(A))+" 0"*(p-1-int(pari.poldegree(A))))); out.flush(); n+=1
    print(q,"%.1fs"%(time.time()-t),flush=True)
print("done %.0fs"%(time.time()-t0))
