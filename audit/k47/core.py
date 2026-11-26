#!/usr/bin/env python3
# Audit k_47 (copied from audit/k41_k43/core.py; added: the report when the gcd is a power of p).
# Audit k_41, k_43, step 3: the rows left by elim.py generate A^n (n = number of columns left).
# Blocks of n rows are chosen greedily (rows ordered by the trace form; a row is taken if it raises the rank of the block modulo the
# prime ideal (l0, x - w0) of A).  The choice only matters for speed.  For each block C: D = Nm(det C) by nmdet2.nm_crt (exact).
# PROOF: the index of the span of all rows divides every nonzero |D|; if the gcd of the |D| is 1, the rows generate A^n.
# Output: out/core_TAG.pkl (blocks and D), out/T_own_P.txt, text on stdout.
#   python3 core.py p RULE NBLOCKS [SECONDS]
import sys, pickle, time, math, os
from nmdet2 import nm_crt, trace_form, zmatrix, hadamard_bound, need, pari
p=int(sys.argv[1]); rule=sys.argv[2]; tag="%d_%s"%(p,rule.replace(":","")); t0=time.time(); sys.set_int_max_str_digits(1000000)
budget=float(sys.argv[4]) if len(sys.argv)>4 else 1e9
st=pickle.load(open("out/elim_%s.pkl"%tag,"rb")); need(st["stuck"] or not st["cols"],"the elimination is not finished")
cols=st["cols"]; n=len(cols); alive=st["alive"]
pad=lambda v:list(v)+[0]*(p-1-len(v))
R={i:[pad(st["rows"][i].get(q,[])) for q in cols] for i in alive}
need(all(set(st["rows"][i])<=set(cols) for i in alive),"support")
print("k_%d, rule %s: core of %d columns, %d rows; pivots before: %d; |T| = %d"%(p,rule,n,len(alive),st["npiv"],len(st["T"])),flush=True)
out=[]; g=0; used=set()
if os.path.exists("out/core_%s.pkl"%tag):
    prev=pickle.load(open("out/core_%s.pkl"%tag,"rb")); out=prev["blocks"]; g=prev["gcd"]; used=set(i for b in out for i in b["rows"])
    print("resuming: %d blocks done, gcd so far %s"%(len(out),"= %d"%g if g<10**9 else "has %d digits"%len(str(g))))
if n>0 and g!=1:
    order=sorted(alive,key=lambda i:(trace_form(R[i],p),i))
    l0=2*p+1
    while not pari.isprime(l0): l0+=2*p
    w0=next(pow(a,(l0-1)//p,l0) for a in range(2,l0) if pow(a,(l0-1)//p,l0)!=1); pw=[pow(w0,d,l0) for d in range(p-1)]
    ev={i:[sum(c*pw[d] for d,c in enumerate(e))%l0 for e in R[i]] for i in alive}
    def greedy(cands):
        ech=[]; piv=[]; chosen=[]
        for i in cands:
            v=ev[i][:]
            for r,c in zip(ech,piv):
                if v[c]:
                    f=v[c]; v=[(a-f*b)%l0 for a,b in zip(v,r)]
            c=next((k for k,a in enumerate(v) if a),None)
            if c is None: continue
            iv=pow(v[c],-1,l0); ech.append([a*iv%l0 for a in v]); piv.append(c); chosen.append(i)
            if len(chosen)==n: break
        return chosen
    while len(out)<int(sys.argv[3]) and time.time()-t0<budget:
        blk=greedy([i for i in order if i not in used])
        if len(blk)<n: print("block %d: only %d rows independent modulo (%d, x - %d) are left"%(len(out)+1,len(blk),l0,w0)); break
        used.update(blk); C=[R[i] for i in blk]
        bound=hadamard_bound(zmatrix(C,p)); D,k=nm_crt(C,p,bound,flint_checks=2); g=math.gcd(g,abs(D)); out.append({"rows":blk,"D":D})
        print("block %d: Nm(det) has %d digits (bound %d digits, %d auxiliary primes); gcd so far %s; %.0f s"%(len(out),len(str(abs(D))),len(str(bound)),k,"= %d"%g if g<10**9 else "has %d digits"%len(str(g)),time.time()-t0),flush=True)
        pickle.dump({"cols":cols,"blocks":out,"gcd":g},open("out/core_%s.pkl"%tag,"wb"))
        if g==1: break
if n==0: g=1
if g==1:
    print("OWN SPAN PASS k_%d (%s): the rows generate A^%d, hence L = A^T for the %d primes of T (smallest %s, largest %d; %d below 10^6)"%(p,rule,n,len(st["T"]),st["T"][:6],st["T"][-1],sum(1 for q in st["T"] if q<10**6)))
    open("out/T_own_%s.txt"%tag,"w").write("".join("%d\n"%q for q in st["T"]))
elif g>1 and p**round(math.log(g,p))==g and g<10**9:
    print("OWN RESULT k_%d (%s): the gcd of the indices is %d = %d^%d.  The rows generate a submodule of A^T whose index divides %d (|T| = %d)."%(p,rule,g,p,round(math.log(g,p)),g,len(st["T"])))
    open("out/T_own_%s.txt"%tag,"w").write("".join("%d\n"%q for q in st["T"]))
else: print("gcd so far %s: not 1.  (Either the run is not finished - run again to continue - or it is a control run RULE+ctlQ, for which a multiple of p = %d is expected.)"%("= %d"%g if g<10**9 else "has %d digits"%len(str(g)),p))
