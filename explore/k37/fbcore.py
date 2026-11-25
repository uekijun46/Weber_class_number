# Exploration (untrusted): reduce the relation matrix over Z[G]/(norm) = Z[zeta_p] by pivots that are units +-sigma^j.
# The module  Z[zeta]^T / (rows)  does not change (up to isomorphism) when a row with a unit at the column q is used to
# clear the column q and then this row and this column are removed.
#   python3 fbcore.py p MAT OUT [stop when this many columns remain]
import sys, pickle, time, numpy as np
p=int(sys.argv[1]); N=p; cols,rows0=pickle.load(open(sys.argv[2],"rb")); stop=int(sys.argv[4]) if len(sys.argv)>4 else 0
def vec(d):
    v=np.zeros(N,dtype=object)
    for j,a in d.items(): v[j%N]+=a
    return v
def red(v): return v-v[N-1]                  # the representative with last coefficient 0
def mono(v):                                 # +-sigma^j modulo the norm element?  returns (sign, j) or None
    w=v-np.sort(v)[N//2]; nz=np.nonzero(w)[0]
    return (int(w[nz[0]]),int(nz[0])) if len(nz)==1 and abs(w[nz[0]])==1 else None
def mul(a,b):
    c=np.zeros(N,dtype=object)
    for j in np.nonzero(a)[0]: c+=a[j]*np.roll(b,j)
    return red(c)
if rows0 and isinstance(rows0[0],dict):
    rows=rows0; keep=int(sys.argv[5]) if len(sys.argv)>5 else len(rows)
    rows=sorted(rows,key=lambda r:sum(int(abs(x).sum()) for x in r.values()))[:keep]
else:
    rows=[{q:red(vec(d)) for q,d in v.items()} for e,v in rows0]
    rows=[{q:x for q,x in r.items() if x.any()} for r in rows]
colrows={}
for i,r in enumerate(rows):
    for q in r: colrows.setdefault(q,set()).add(i)
alive=set(range(len(rows))); t0=time.time(); npiv=0; dropped=[]
def drop_row(i):
    for q in rows[i]:
        colrows[q].discard(i)
    alive.discard(i)
def prune():                                 # primes that occur in one row only: remove the prime and the row
    ch=True
    while ch:
        ch=False
        for q in [q for q,s in colrows.items() if len(s)<=1]:
            for i in list(colrows[q]): drop_row(i)
            del colrows[q]; dropped.append(q); ch=True
prune()
print("start: %d rows, %d columns"%(len(alive),len(colrows)),flush=True)
while len(colrows)>stop:
    best=None
    for q,s in sorted(colrows.items(),key=lambda t:len(t[1]))[:40]:
        for i in s:
            m=mono(rows[i][q])
            if m:
                cost=(len(rows[i])-1)*(len(s)-1)
                if best is None or cost<best[0]: best=(cost,q,i,m)
    if best is None:
        for q,s in colrows.items():
            for i in s:
                m=mono(rows[i][q])
                if m:
                    cost=(len(rows[i])-1)*(len(s)-1)
                    if best is None or cost<best[0]: best=(cost,q,i,m)
    if best is None: break
    cost,q,i,(sg,j)=best; r=rows[i]; inv=np.zeros(N,dtype=object); inv[(-j)%N]=sg; inv=red(inv)
    for k in list(colrows[q]):
        if k==i: continue
        f=mul(rows[k][q],inv); rk=rows[k]
        for q2,x in r.items():
            if q2==q: continue
            y=rk.get(q2); z=(y if y is not None else 0)-mul(f,x)
            if z.any():
                rk[q2]=z; colrows[q2].add(k)
            elif y is not None:
                del rk[q2]; colrows[q2].discard(k)
        del rk[q]
    drop_row(i); del colrows[q]; npiv+=1
    for k in [k for k in alive if not rows[k]]: alive.discard(k)
    if npiv%(100 if len(colrows)>100 else 10)==0:
        print("pivots %d: %d rows, %d columns, entries %d, largest coefficient %d (%.0fs)"%(npiv,len(alive),len(colrows),sum(len(rows[k]) for k in alive),max(int(abs(x).max()) for k in alive for x in rows[k].values()),time.time()-t0),flush=True)
    if time.time()-t0>500: break
T=sorted(colrows)
print("end: pivots %d; core %d rows, %d columns (smallest primes %s); entries %d; unit entries left %d; primes removed as single: %d (%.0fs)"%(npiv,len(alive),len(T),T[:8],sum(len(rows[k]) for k in alive),sum(1 for k in alive for x in rows[k].values() if mono(x)),len(dropped),time.time()-t0))
pickle.dump((T,[rows[k] for k in sorted(alive)]),open(sys.argv[3],"wb"))
