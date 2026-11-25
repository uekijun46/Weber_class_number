#!/usr/bin/env python3
# Audit k_37: the audit's own proof that relation rows generate A^T, A = Z[x]/(Phi_p), as functions (the same algorithm as
# elim_own.py + core_own.py, for any odd prime p).  Nothing from relations/ is imported.
#   rows: list of dicts {q: [p-1 integers]} (entries reduced modulo Phi_p, coefficients of 1, x, ..., x^(p-2)); T = all q that occur.
#   eliminate: pivots on entries +-x^j only, smallest Markowitz cost; no row is dropped, no prime is removed.
#   core_gcd:  blocks of n rows of the remaining n columns; D = Nm(det) by nmdet.nm_crt; returns the gcd of the |D| (0 if no block has
#              nonzero determinant).  n = 0 gives 1.
# own_span returns (T, gcd).  gcd = 1  <=>  proved that the rows generate A^T.  In general the index of the span is a divisor of gcd... no:
# the index [A^T : L] DIVIDES every nonzero |D|, hence divides the gcd; if a prime ideal of T is not principal, the gcd cannot be 1.
import sys, math, time
from flint import fmpz_poly
from nmdet import nm_crt, trace_form, zmatrix, hadamard_bound, need, pari
def eliminate(rows_in,p,verbose=False):
    Phi=fmpz_poly([1]*p); ZERO=fmpz_poly([]); ONE=fmpz_poly([1]); red=lambda a:a%Phi
    def unit(a):
        co=[int(v) for v in a.coeffs()]; nz=[j for j,v in enumerate(co) if v]
        if len(nz)==1 and abs(co[nz[0]])==1: return co[nz[0]],nz[0]
        if p>2 and len(co)==p-1 and len(set(co))==1 and abs(co[0])==1: return -co[0],p-1
        return None
    rows=[{q:red(fmpz_poly(v)) for q,v in r.items()} for r in rows_in]
    rows=[{q:v for q,v in r.items() if v!=ZERO} for r in rows]
    alive=set(i for i,r in enumerate(rows) if r); T=sorted(set(q for r in rows for q in r)); cols=set(T); col={q:set() for q in cols}; npiv=0; t0=time.time()
    for i in alive:
        for q in rows[i]: col[q].add(i)
    while cols:
        best=None; found=0
        for q in sorted(cols,key=lambda q:(len(col[q]),q)):
            if found>=60: break
            hit=False
            for i in col[q]:
                u=unit(rows[i][q])
                if u is not None:
                    hit=True; cost=((len(rows[i])-1)*(len(col[q])-1),q,i)
                    if best is None or cost<best[0]: best=(cost,u)
            if hit: found+=1
            if best is not None and best[0][0]==0: break
        if best is None: break
        (cost,q,i),(sg,j)=best; r=rows[i]; inv=red(fmpz_poly([0]*((p-j)%p)+[sg]))
        need(red(r[q]*inv)==ONE,"pivot inverse")
        others=[(q2,x) for q2,x in r.items() if q2!=q]
        for k in list(col[q]):
            if k==i: continue
            rk=rows[k]; fct=red(rk[q]*inv)
            for q2,x in others:
                y=rk.get(q2,ZERO)-red(fct*x)
                if y!=ZERO: rk[q2]=y; col[q2].add(k)
                elif q2 in rk: del rk[q2]; col[q2].discard(k)
            del rk[q]
            if not rk: alive.discard(k)
        for q2,x in others: col[q2].discard(i)
        alive.discard(i); del col[q]; cols.discard(q); npiv+=1
        if verbose and npiv%500==0: print("  pivots %d; columns %d; rows %d; %.0f s"%(npiv,len(cols),len(alive),time.time()-t0),flush=True)
    need(all(set(rows[i])<=cols for i in alive),"support")
    cl=sorted(cols); pad=lambda v:[int(c) for c in v.coeffs()]+[0]*(p-1-len(v.coeffs()))
    left={i:[pad(rows[i][q]) if q in rows[i] else [0]*(p-1) for q in cl] for i in sorted(alive)}
    return T,cl,left,npiv
def core_gcd(cl,left,p,maxblocks=6,verbose=False):
    n=len(cl)
    if n==0: return 1,[]
    order=sorted(left,key=lambda i:(trace_form(left[i],p),i)); g=0; info=[]; used=set(); t0=time.time()
    l0=100; tries=0
    while len(info)<maxblocks and tries<4:
        l0+=1
        while l0%p!=1 or not pari.isprime(l0): l0+=1
        tries+=1
        w0=next(pow(a,(l0-1)//p,l0) for a in range(2,l0) if pow(a,(l0-1)//p,l0)!=1); pw=[pow(w0,d,l0) for d in range(p-1)]
        ev={i:[sum(c*pw[d] for d,c in enumerate(e))%l0 for e in left[i]] for i in left}
        while len(info)<maxblocks:
            ech=[]; piv=[]; blk=[]
            for i in order:
                if i in used: continue
                v=ev[i][:]
                for r,c in zip(ech,piv):
                    if v[c]:
                        f=v[c]; v=[(a-f*b)%l0 for a,b in zip(v,r)]
                c=next((k for k,a in enumerate(v) if a),None)
                if c is None: continue
                iv=pow(v[c],-1,l0); ech.append([a*iv%l0 for a in v]); piv.append(c); blk.append(i)
                if len(blk)==n: break
            if len(blk)<n: break                       # no further block that is independent modulo (l0, x - w0); try another l0
            used.update(blk); C=[left[i] for i in blk]; bound=hadamard_bound(zmatrix(C,p)); D,k=nm_crt(C,p,bound); g=math.gcd(g,abs(D))
            info.append((abs(D).bit_length(),bound.bit_length(),k,g if g<10**12 else -1))
            if verbose: print("  block %d: Nm(det) has %d digits (bound %d digits, %d auxiliary primes); gcd so far %s; %.0f s"%(len(info),len(str(abs(D))),len(str(bound)),k,"= %d"%g if g<10**12 else "has %d digits"%len(str(g)),time.time()-t0),flush=True)
            if g==1: return 1,info
    return g,info
def own_span(rows,p,maxblocks=6,verbose=False):
    T,cl,left,npiv=eliminate(rows,p,verbose)
    if verbose: print("  elimination: %d pivots; %d columns and %d rows left"%(npiv,len(cl),len(left)),flush=True)
    g,info=core_gcd(cl,left,p,maxblocks,verbose)
    return T,g,(npiv,len(cl),len(left),info)
if __name__=="__main__":                              # k_37:  python3 ownspan.py out/rows.json
    import json
    sys.set_int_max_str_digits(1000000)
    p=37; raw=json.load(open(sys.argv[1])); rows=[]
    for r in raw:
        d={}
        for q,J in r:
            e=[0]*p
            for j in J: e[j]=1
            d[q]=[e[k]-e[p-1] for k in range(p-1)]
        rows.append(d)
    T,g,st=own_span(rows,p,6,True)
    need(g==1,"gcd = %d"%g)
    print("OWN SPAN PASS (ownspan.py): L = A^T for the %d primes of T (smallest %s, largest %d); pivots %d, core %d columns"%(len(T),T[:6],T[-1],st[0],st[1]))
