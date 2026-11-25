#!/usr/bin/env python3
# Audit k_37: for the matrices C recorded from span.py (out/span_blocks.pkl), the determinant of the integer matrix zmatrix(C)
# (size 36 n = 4320) modulo primes near 2^63 that are not 1 mod 37, compared with +-D of span.py.  A consistency check of
# "index = |Nm(det C)|" by integer linear algebra (flint nmod_mat.det), not a proof.     python3 zdet_span.py BLOCK NPRIMES
import sys, pickle, time
from flint import nmod_mat, fmpz
from nmdet import zmatrix, need
p=37; B=pickle.load(open("out/span_blocks.pkl","rb"))[int(sys.argv[1])]; n=B["n"]; fl=B["flat"]; D=B["D"]
C=[[fl[(i*n+j)*(p-1):(i*n+j+1)*(p-1)] for j in range(n)] for i in range(n)]
Z=zmatrix(C,p); m=len(Z); flat=[c for r in Z for c in r]; t0=time.time(); l=2**63
for _ in range(int(sys.argv[2])):
    l-=1
    while not fmpz(l).is_prime() or l%p==1: l-=1
    dz=int(nmod_mat(m,m,flat,l).det()); need(dz in (D%l,(-D)%l),"mismatch modulo %d"%l)
    print("set %d of span.py: det of the %d x %d integer matrix = %sD modulo %d; %.0f s"%(int(sys.argv[1])+1,m,m,"" if dz==D%l else "-",l,time.time()-t0),flush=True)
