#!/usr/bin/env python3
# Audit k_41, k_43: supporting evidence (and, for each single prime ideal, an exact statement): ranks modulo prime ideals of A.
#   For a prime l = 1 mod p and w of order p modulo l, (l, x - w^k), k = 1..p-1, are the prime ideals of A above l (residue field F_l);
#   (p, x - 1) is the prime ideal above p (residue field F_p).  If the matrix of the relation rows has full column rank modulo a prime
#   ideal P, then P does not divide the index ideal of L in A^T (Nakayama), i.e. the quotient A^T/L has no quotient A/P.
#   mode "core": the rows left by elim.py (n columns);  mode "full": the rows of out/rows_P.json before any elimination (|T| columns).
#   python3 rank_check.py p core RULE LMAX      |      python3 rank_check.py p full l [k1 k2 ...]
import sys, json, pickle, time
from flint import nmod_mat
from fieldp import need, is_prime_td
p=int(sys.argv[1]); mode=sys.argv[2]; t0=time.time()
def root(l):
    w=next(pow(a,(l-1)//p,l) for a in range(2,l) if pow(a,(l-1)//p,l)!=1); need(pow(w,p,l)==1 and w!=1,"w"); return w
if mode=="core":
    rule=sys.argv[3]; st=pickle.load(open("out/elim_%d_%s.pkl"%(p,rule),"rb")); cols=st["cols"]; n=len(cols); alive=st["alive"]
    ent=[[st["rows"][i].get(q,[]) for q in cols] for i in alive]
    pad=lambda v:list(v)+[0]*(p-1-len(v)); flat=[c for r in ent for e in r for c in pad(e)]; nr=len(ent)
    r=nmod_mat(nr,n,[sum(e)%p for rw in ent for e in rw],p).rank()        # x = 1 modulo p
    print("k_%d core (%s): %d columns, %d rows; rank modulo (p, x - 1): %d"%(p,rule,n,nr,r)); bad=[] if r==n else [(p,0)]
    ls=[l for l in range(p+1,int(sys.argv[4])) if l%p==1 and is_prime_td(l)]
    for l in ls:
        w=root(l); E=nmod_mat(nr*n,p-1,flat,l); V=nmod_mat(p-1,p-1,[pow(w,(k+1)*d,l) for d in range(p-1) for k in range(p-1)],l)
        EV=(E*V).entries(); rk=[nmod_mat(nr,n,EV[k::p-1],l).rank() for k in range(p-1)]      # column k: evaluation at w^(k+1)
        if min(rk)<n: bad+=[(l,k+1) for k,v in enumerate(rk) if v<n]
    print("primes l = 1 mod %d below %s: %s"%(p,sys.argv[4],ls))
    print("rank modulo (l, x - w^k) is %d for all these l and all k = 1..%d: %s%s; %.0f s"%(n,p-1,not bad,"" if not bad else "; exceptions: %s"%bad,time.time()-t0))
else:
    l=int(sys.argv[3]); ks=[int(a) for a in sys.argv[4:]] or list(range(1,p)); raw=json.load(open("out/rows_%d.json"%p))
    T=sorted(set(q for r in raw for q,J in r)); idx={q:t for t,q in enumerate(T)}; res=[]
    if l==p: pts=[1]
    else: w=root(l); pts=[pow(w,k,l) for k in ks]
    for pt in pts:
        pw=[pow(pt,j,l) for j in range(p)]; Mx=nmod_mat(len(raw),len(T),l)
        for i,r in enumerate(raw):
            for q,J in r: Mx[i,idx[q]]=sum(pw[j] for j in J)%l
        res.append(Mx.rank()); print("  l = %d, x = %d: rank %d of %d; %.0f s"%(l,pt,res[-1],len(T),time.time()-t0),flush=True)
    print("k_%d full matrix (%d rows, %d columns) modulo l = %d at %d points: all ranks equal to %d: %s"%(p,len(raw),len(T),l,len(pts),len(T),all(v==len(T) for v in res)))
