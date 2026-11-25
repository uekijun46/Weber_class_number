#!/usr/bin/env python3
# Audit k_37: for the blocks of out/core_blocks.pkl, compare D = Nm(det C) (nmdet.nm_crt) with the determinant of the integer matrix
# zmatrix(C) of size 36 n.   Mode "mod K": det modulo K primes near 2^63 that are NOT 1 mod 37 (flint nmod_mat.det): a consistency check.
# Mode "exact": flint fmpz_mat.det (exact).   python3 zdet_check.py BLOCK mod K | exact
import sys, pickle, time
from flint import fmpz_mat, nmod_mat, fmpz
from nmdet import zmatrix, need
sys.set_int_max_str_digits(1000000)
p=37; d=pickle.load(open(sys.argv[4] if len(sys.argv)>4 else "out/core_blocks.pkl","rb")); b=d["blocks"][int(sys.argv[1])]; C=b["C"]; D=b["D"]
Z=zmatrix(C,p); m=len(Z); flat=[c for r in Z for c in r]; t0=time.time()
if sys.argv[2]=="exact":
    dz=int(fmpz_mat(m,m,flat).det()); print("block %s: exact determinant of the %d x %d integer matrix: %d digits; equal to Nm(det C): %s; %.0f s"%(sys.argv[1],m,m,len(str(abs(dz))),dz==D,time.time()-t0)); need(abs(dz)==abs(D),"mismatch")
else:
    l=2**63; ok=0
    for _ in range(int(sys.argv[3])):
        l-=1
        while not fmpz(l).is_prime() or l%p==1: l-=1
        dz=int(nmod_mat(m,m,flat,l).det()); need(dz==D%l,"mismatch modulo %d"%l); ok+=1
        print("  block %s: det of the %d x %d integer matrix = Nm(det C) modulo %d; %.0f s"%(sys.argv[1],m,m,l,time.time()-t0),flush=True)
    print("block %s: agreement modulo %d primes"%(sys.argv[1],ok))
