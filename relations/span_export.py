# Export (untrusted search data -> certificate for span.py).  Input: output files of span_search.py.
# Keeps the relations all of whose primes occur in at least two kept relations, and drops relations in which a prime ideal
# occurs with exponent > 1.   python3 span_export.py p OUT FILE [FILE ...]
import sys, pickle, random
from field import Field
p,out=int(sys.argv[1]),sys.argv[2]; F=Field(p,1); N=F.N; rnd=random.Random(1)
rels=[]
for fn in sys.argv[3:]: rels+=pickle.load(open(fn,"rb"))
def factor(x,lp,primes):
    fac={}
    for r in ([p]+([lp] if lp else [])):
        while x%r==0: x//=r; fac[r]=fac.get(r,0)+1
    for q in primes:
        if x==1: break
        while x%q==0: x//=q; fac[q]=fac.get(q,0)+1
    return fac if x==1 else None
import flint
Y=int(float(sys.argv[sys.argv.index("--Y")+1])) if "--Y" in sys.argv else 10**6
fb=[q for q in range(2,Y) if F.split(q) and flint.fmpz(q).is_prime()]
cnt={}
for r in rels:
    if r[3]: cnt[r[3]]=cnt.get(r[3],0)+1
Rq={}; rows=[]
for c,coef,x,lp in rels:
    if lp and cnt[lp]<2: continue
    fac=factor(x,lp,fb)
    if fac is None: continue
    ok=True
    for q,k in fac.items():
        if q==p: continue
        if q not in Rq: Rq[q]=F.periods_mod(q,rnd)
        R=Rq[q]
        if sum(1 for j in range(N) if (c+sum(a*R[(i+j)%N] for i,a in coef))%q==0)!=k: ok=False; break
    if ok: rows.append((c,coef,fac))
alive=set(range(len(rows))); ch=True
while ch:
    occ={}
    for i in alive:
        for q in rows[i][2]:
            if q!=p: occ[q]=occ.get(q,0)+1
    bad={i for i in alive if any(q!=p and occ[q]<2 for q in rows[i][2])}; ch=bool(bad); alive-=bad
with open(out,"w") as f:
    f.write("# field %d 1\n# relations for span.py: c | i:a ... | factorization of |N(c + sum a_i eta_i)|\n"%p)
    for i in sorted(alive):
        c,coef,fac=rows[i]
        f.write("%d | %s | %s\n"%(c," ".join("%d:%d"%(i_,a) for i_,a in coef)," ".join("%d^%d"%(q,k) for q,k in sorted(fac.items()))))
print("relations read %d, written %d"%(len(rels),len(alive)))
