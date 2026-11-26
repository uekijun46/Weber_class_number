#!/usr/bin/env python3
# Audit k_47: where does the factor 47 come from?  No elimination over A is used here.
# The only prime ideal of A = Z[x]/(Phi_47) above 47 is (1 - x), and A/(1 - x) = F_47 (x -> 1).  Under x -> 1 the entry sum_{j in J} x^j
# becomes |J| (the exponent of q in the norm).  So the matrix of all rows modulo (1 - x) is the integer matrix of the exponents modulo 47.
# A^T/L has a quotient of order 47 iff this matrix has rank < |T| over F_47; (A^T/L)/(1 - x) = F_47^T / (row space), of dimension |T| - rank.
# The right kernel {v : M v = 0} is the space of F_47-linear functionals on F_47^T that vanish on the rows.
#   python3 rank47.py [dropA,B]         (data lines to leave out; their primes stay in T unless they occur nowhere else)
import sys, json, time
from flint import nmod_mat
p=47; t0=time.time(); raw=json.load(open("out/rows_47.json"))
drop=[int(v) for v in sys.argv[1][4:].split(",")] if len(sys.argv)>1 else []
for i in drop: raw[i]=[]
T=sorted(set(q for r in raw for q,J in r)); idx={q:t for t,q in enumerate(T)}
M=nmod_mat(len(raw),len(T),p)
for i,r in enumerate(raw):
    for q,J in r: M[i,idx[q]]=len(J)%p
print("matrix of exponents modulo 47: %d rows (%d left out), %d columns"%(len(raw),len(drop),len(T)),flush=True)
X,nul=M.nullspace()
print("rank over F_47: %d; |T| - rank = %d; %.0f s"%(len(T)-nul,nul,time.time()-t0),flush=True)
for b in range(nul):
    v=[int(X[t,b]) for t in range(len(T))]
    # check M v = 0 directly from the rows (independent of the nullspace routine)
    ok=all(sum(len(J)*v[idx[q]] for q,J in r)%p==0 for r in raw)
    sup=[(T[t],v[t]) for t in range(len(T)) if v[t]]
    print("kernel vector %d: M v = 0 checked row by row: %s; support (prime, coefficient): %s"%(b,ok,sup if len(sup)<=20 else "%d primes"%len(sup)))
