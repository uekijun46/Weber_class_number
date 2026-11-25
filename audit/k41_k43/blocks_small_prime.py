#!/usr/bin/env python3
# Audit k_41: why the gcd of span.py passes through 83.  For each set C recorded from span.py (span_rerun.py): the exponent of l in
# D = |Nm(det C)| and the list of k for which det C(w^k) = 0 modulo l (w of order p modulo l), i.e. the prime ideals (l, x - w^k) of A
# that divide det C.  If different sets are divisible by different prime ideals above l, the common factor l of the numbers D is not a
# common factor of the ideals (det C), and nothing is wrong.     python3 blocks_small_prime.py p l
import sys, pickle
from flint import nmod_mat
sys.set_int_max_str_digits(1000000)
p=int(sys.argv[1]); l=int(sys.argv[2]); blocks=pickle.load(open("out/span_blocks_%d.pkl"%p,"rb"))
w=next(pow(a,(l-1)//p,l) for a in range(2,l) if pow(a,(l-1)//p,l)!=1)
for b,B in enumerate(blocks):
    n=B["n"]; fl=B["flat"]; D=B["D"]; v=0
    while D and D%l==0: D//=l; v+=1
    ks=[]
    for k in range(1,p):
        pw=[pow(w,k*d,l) for d in range(p-1)]
        Mx=nmod_mat(n,n,[sum(c*x for c,x in zip(fl[e*(p-1):(e+1)*(p-1)],pw))%l for e in range(n*n)],l)
        if int(Mx.det())==0: ks.append((k,n-Mx.rank()))
    print("set %d: exponent of %d in |Nm(det)|: %d; (k, corank of C(w^k) mod %d) for the k with det = 0: %s"%(b+1,l,v,l,ks))
