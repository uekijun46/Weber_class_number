# Exploration (untrusted): from the relations of fbfast.py to a matrix over Z[G], G = Gal(k_p/Q) = <tau>, and rank tests.
# For a split prime q fix phi_0: O -> F_q (eta_i -> R[i]) and put q_j = ker(phi_0 tau^j) = sigma^j(q_0), sigma = tau^{-1}.
# If q^k || N(alpha) and exactly k of the phi_j(alpha) vanish, j in J, then the q-part of (alpha) is (sum_{j in J} sigma^j) q_0.
# Rows with the same large prime are subtracted from each other after multiplication by a power of sigma.
#   python3 fbmat.py p Y OUT FILE [FILE ...]
import sys, pickle, random, time, flint
sys.path.insert(0,"../../relations")
from field import Field
p,Y,out=int(sys.argv[1]),int(float(sys.argv[2])),sys.argv[3]; F=Field(p,1); N=F.N
rels=[]
for fn in sys.argv[4:]: rels+=pickle.load(open(fn,"rb"))
fb=[q for q in range(2,Y) if F.split(q) and flint.fmpz(q).is_prime()]
rnd=random.Random(1); Rq={}
def per(q):
    if q not in Rq: Rq[q]=F.periods_mod(q,rnd)
    return Rq[q]
def pattern(c,coef,q):
    R=per(q); return [j for j in range(N) if (c+sum(a*R[(i+j)%N] for i,a in coef))%q==0]
cnt={}
for r in rels:
    if r[3]: cnt[r[3]]=cnt.get(r[3],0)+1
t0=time.time(); full=[]; part={}; bad=0
for c,coef,x,lp in rels:
    if lp and cnt[lp]<2: continue
    v={}; ok=True; y=x
    while y%p==0: y//=p
    if lp: y//=lp
    for q in fb:
        if y%q==0:
            k=0
            while y%q==0: y//=q; k+=1
            J=pattern(c,coef,q)
            if len(J)!=k: ok=False; break
            v[q]={j:1 for j in J}
            if y==1: break
    if not ok or y!=1: bad+=1; continue
    if lp is None: full.append(((c,coef,0),v))
    else:
        J=pattern(c,coef,lp)
        if len(J)!=1: bad+=1; continue
        s=J[0]; part.setdefault(lp,[]).append(((c,coef,s),{q:{(j-s)%N:a for j,a in d.items()} for q,d in v.items()}))
rows=[(( 1,e),) +(v,) for e,v in full]      # a row: (list of (sign, element with shift)), exponents
rows=[([(1,e)],v) for e,v in full]
for lp,vs in part.items():
    e0,v0=vs[0]
    for e,v in vs[1:]:
        d={q:dict(x) for q,x in v.items()}
        for q,x in v0.items():
            dq=d.setdefault(q,{})
            for j,a in x.items():
                dq[j]=dq.get(j,0)-a
                if dq[j]==0: del dq[j]
        d={q:x for q,x in d.items() if x}
        if d: rows.append(([(1,e),(-1,e0)],d))
print("relations %d; discarded %d; full %d; large primes used %d; rows %d (%.0fs)"%(len(rels),bad,len(full),len(part),len(rows),time.time()-t0),flush=True)
def canon(v):
    q0=min(v); s=min(v[q0]); sg=1 if v[q0][s]>0 else -1
    return tuple(sorted((q,tuple(sorted(((j-s)%N,a*sg) for j,a in d.items()))) for q,d in v.items()))
uniq={}
for e,v in rows: uniq.setdefault(canon(v),(e,v))
rows=list(uniq.values()); cols=sorted({q for e,v in rows for q in v}); pos={q:i for i,q in enumerate(cols)}
print("distinct rows %d; primes of the factor base that occur: %d of %d; missing below 6000: %s"%(len(rows),len(cols),len(fb),[q for q in fb if q<6000 and q not in pos]),flush=True)
pickle.dump((cols,rows),open(out,"wb"))
for l in [q for q in range(2,20000) if q%N==1 and flint.fmpz(q).is_prime()][:2]:
    w0=next(pow(x,(l-1)//N,l) for x in range(2,l) if pow(x,(l-1)//N,l)!=1)
    for k in (1,2):
        w=pow(w0,k,l); pw=[pow(w,j,l) for j in range(N)]; M=flint.nmod_mat(len(rows),len(cols),l)
        for i,(e,v) in enumerate(rows):
            for q,d in v.items(): M[i,pos[q]]=sum(a*pw[j] for j,a in d.items())%l
        t=time.time(); print("l=%d, embedding %d: rank %d of %d (%.0fs)"%(l,k,M.rank(),len(cols),time.time()-t),flush=True)
