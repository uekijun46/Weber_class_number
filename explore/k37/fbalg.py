# Exploration (untrusted): do the relations generate everything?  Rank of the relation matrix over F_l for embeddings zeta -> w.
# Rows are elements of Z[G]^T (G = Gal(k_p/Q), cyclic of order p, T = set of split primes); modulo the norm element the ring is
# Z[zeta_p].  If the rows generate Z[zeta_p]^T, every prime ideal above a prime of T is principal.
#   python3 fbalg.py p FILE [FILE ...]
import sys, pickle, random, flint
sys.path.insert(0,"../../relations")
from field import Field
p=int(sys.argv[1]); F=Field(p,1); N=F.N
rels=[]
for fn in sys.argv[2:]: rels+=pickle.load(open(fn,"rb"))
random.seed(1); Rq={}
def per(q):
    if q not in Rq: Rq[q]=F.periods_mod(q)
    return Rq[q]
def pattern(c,coef,q):                       # the set of j with alpha in tau^j(q_0)
    R=per(q); return [j for j in range(N) if (c+sum(a*R[(i-j)%N] for i,a in coef))%q==0]
rows=[]; partial={}; bad=0
for c,coef,e,lp in rels:
    v={}; ok=True
    for q,k in e.items():
        J=pattern(c,coef,q)
        if len(J)!=k: ok=False; break
        v[q]={j:1 for j in J}
    if not ok: bad+=1; continue
    if lp is None: rows.append(v)
    else:
        J=pattern(c,coef,lp)
        if len(J)!=1: bad+=1; continue
        s=J[0]; v={q:{(j-s)%N:a for j,a in d.items()} for q,d in v.items()}     # shifted so that the element lies in Q_0
        partial.setdefault(lp,[]).append(v)
comb=0
for lp,vs in partial.items():
    for v in vs[1:]:
        d={q:dict(x) for q,x in vs[0].items()}
        for q,x in v.items():
            dq=d.setdefault(q,{})
            for j,a in x.items(): dq[j]=dq.get(j,0)-a
        rows.append(d); comb+=1
print("relations read: %d; discarded (a prime ideal with valuation > 1): %d; full: %d; from pairs with the same large prime: %d"%(len(rels),bad,len(rows)-comb,comb))
def canon(v):                                # rows that differ by a power of tau or by the sign are the same relation
    q0=min(v); s=min(v[q0]); sg=1 if v[q0][s]>0 else -1
    return tuple(sorted((q,tuple(sorted(((j-s)%N,a*sg) for j,a in d.items() if a))) for q,d in v.items()))
uniq={}
for v in rows:
    v={q:{j:a for j,a in d.items() if a} for q,d in v.items()}; v={q:d for q,d in v.items() if d}
    if v: uniq[canon(v)]=v
print("distinct relations up to conjugation and sign: %d of %d"%(len(uniq),len(rows))); rows=list(uniq.values())
while True:                                  # remove primes that occur in at most one relation, with that relation
    cnt={}
    for v in rows:
        for q in v: cnt[q]=cnt.get(q,0)+1
    single={q for q,k in cnt.items() if k<2}
    if not single: break
    rows=[v for v in rows if not (set(v)&single)]
T=sorted(cnt); pos={q:i for i,q in enumerate(T)}
print("after removing primes that occur once: %d relations, %d primes (smallest %s, largest %s)"%(len(rows),len(T),T[:5],T[-1] if T else None))
if rows and len(rows)>=len(T):
    for l in [q for q in range(2,20000) if q%N==1 and flint.fmpz(q).is_prime()][:2]:
        w0=next(pow(x,(l-1)//N,l) for x in range(2,l) if pow(x,(l-1)//N,l)!=1); ranks=[]
        for k in (1,2,N-1):
            w=pow(w0,k,l); M=flint.nmod_mat(len(rows),len(T),l)
            for i,v in enumerate(rows):
                for q,d in v.items(): M[i,pos[q]]=sum(a*pow(w,j,l) for j,a in d.items())%l
            ranks.append(M.rank())
        print("l=%d: ranks for three embeddings: %s of %d"%(l,ranks,len(T)))
    pickle.dump((T,rows),open("matrix_%d.pkl"%p,"wb"))
