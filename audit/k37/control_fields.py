#!/usr/bin/env python3
# Audit k_37, control experiments: cyclic fields of prime degree d and prime conductor l with their class groups (PARI bnfinit,
# certified by bnfcertify, so the class groups do not depend on GRH).
import sys, cypari2
pari=cypari2.Pari(); pari.allocatemem(10**9,silent=True)
for d,lmax in ((3,1300),(5,2200),(7,2500),(11,1500),(13,1200)):
    out=[]
    for l in range(2,lmax):
        if (l-1)%d or not pari.isprime(l): continue
        bnf=pari.bnfinit(pari("polsubcyclo(%d,%d)"%(l,d)),1)
        if int(bnf[7][0][0])>1 or l<80:
            need_cert=int(pari.bnfcertify(bnf)); out.append((l,int(bnf[7][0][0]),[int(c) for c in bnf[7][0][1]],need_cert))
    print("d=%d:"%d," ".join("%d:h=%d%s%s"%(l,h,cyc,"" if ok==1 else "(not certified)") for l,h,cyc,ok in out),flush=True)
