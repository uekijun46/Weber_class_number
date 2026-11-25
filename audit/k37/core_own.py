#!/usr/bin/env python3
# Audit k_37, step 3: the rows left by elim_own.py generate A^n (n = number of columns left).
# Blocks of n rows are chosen greedily (rows ordered by the trace form, a row is taken if it raises the rank of the block modulo the
# prime ideal (l0, x - w0) of A); the choice only matters for speed.  For each block C: D = Nm(det C) by nmdet.nm_crt (exact).
# PROOF: if the gcd of the |D| is 1, the rows generate A^n.   Output: out/core_blocks.pkl (blocks and D), text on stdout.
#   python3 core_own.py NBLOCKS
import sys, pickle, time, math
from nmdet import nm_crt, trace_form, zmatrix, hadamard_bound, need
p=37; t0=time.time(); sys.set_int_max_str_digits(1000000)
st=pickle.load(open("out/elim_state.pkl","rb")); need(st["stuck"] or not st["cols"],"the elimination is not finished")
cols=st["cols"]; n=len(cols); alive=st["alive"]
pad=lambda v:list(v)+[0]*(p-1-len(v))
R={i:[pad(st["rows"][i].get(q,[])) for q in cols] for i in alive}
need(all(set(st["rows"][i])<=set(cols) for i in alive),"support")
print("core: %d columns, %d rows; pivots before: %d; |T| = %d"%(n,len(alive),st["npiv"],len(st["T"])),flush=True)
order=sorted(alive,key=lambda i:(trace_form(R[i],p),i))
l0=149; w0=next(pow(a,(l0-1)//p,l0) for a in range(2,l0) if pow(a,(l0-1)//p,l0)!=1); pw=[pow(w0,d,l0) for d in range(p-1)]
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
used=set(); g=0; out=[]
for b in range(int(sys.argv[1])):
    blk=greedy([i for i in order if i not in used])
    if len(blk)<n: print("block %d: only %d independent rows left"%(b+1,len(blk))); break
    used.update(blk); C=[R[i] for i in blk]
    bound=hadamard_bound(zmatrix(C,p)); D,k=nm_crt(C,p,bound); g=math.gcd(g,abs(D)); out.append({"rows":blk,"C":C,"D":D})
    print("block %d: Nm(det) has %d digits (bound %d digits, %d auxiliary primes); gcd so far %s; %.0f s"%(b+1,len(str(abs(D))),len(str(bound)),k,"= %d"%g if g<10**9 else "has %d digits"%len(str(g)),time.time()-t0),flush=True)
    pickle.dump({"cols":cols,"blocks":out,"gcd":g},open("out/core_blocks.pkl","wb"))
    if g==1: break
need(g==1,"gcd is not 1")
print("OWN SPAN PASS: the rows generate A^%d, hence (elim_own.py) L = A^T for the %d primes of T (smallest %s, largest %d)"%(n,len(st["T"]),st["T"][:6],st["T"][-1]))
open("out/T_own.txt","w").write("".join("%d\n"%q for q in st["T"]))
