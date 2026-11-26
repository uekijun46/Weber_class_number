#!/usr/bin/env python3
# Audit k_47: does the changed stopping rule of relations/span.py change the recorded runs for k_37, k_41, k_43?
# Input: the numbers D = |Nm(det C)| of step (c) recorded from span.py by the earlier audits
#   (audit/k37/out/span_blocks.pkl, audit/k41_k43/out/compare_blocks_41.json, compare_blocks_43.json).
# For every set: old rule "stop if g == 1", new rule "stop if g != 0 and g is a power of p".  They must agree at every step.
import json, pickle, math, sys, os
sys.set_int_max_str_digits(0)
here=os.path.dirname(os.path.abspath(__file__))
def ppower(g,p):
    while g%p==0: g//=p
    return g==1
def trace(p,Ds):
    g=0; same=True; out=[]
    for D in Ds:
        g=math.gcd(g,abs(int(D))); old=(g==1); new=(g!=0 and ppower(g,p)); same=same and old==new
        vp=0
        while g and g%p**(vp+1)==0: vp+=1
        out.append("gcd %s (divisible by p exactly %d times), old rule stops: %s, new rule stops: %s"%("= %d"%g if g<10**6 else "has %d digits"%len(str(g)),vp,old,new))
        if old or new: break
    print("k_%d: %d sets"%(p,len(out)))
    for t in out: print("   "+t)
    print("   the two rules agree at every step: %s; final gcd %d"%(same,g))
b=pickle.load(open(os.path.join(here,"..","k37","out","span_blocks.pkl"),"rb")); trace(37,[x["D"] for x in b])
for p in (41,43):
    d=json.load(open(os.path.join(here,"..","k41_k43","out","compare_blocks_%d.json"%p))); trace(p,[d[str(i)] for i in range(len(d))])
