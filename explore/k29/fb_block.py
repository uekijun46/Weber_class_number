# Exploration (untrusted): can the primes above the small split primes of k_29 be shown principal together?
# Factor base FB = split primes < Y.  Elements whose norm is 29^a * prod_{q in FB} q^(k_q) and which lie in exactly k_q distinct
# primes above q (valuations 0 or 1) give 0/1 vectors in Z^(29*|FB|); if the vectors of all conjugates span the whole lattice,
# then every prime above a prime of FB is principal.
import sys, pickle, flint
sys.path.insert(0,"../../relations"); sys.path.insert(0,".")
from field import Field
p,n,Y=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
F=Field(p,n); N=F.N
rows=pickle.load(open("rows_%d_%d.pkl"%(p,n),"rb"))
seen=set(); smooth=[]
for c,coef,fac in rows:
    key=(c,tuple(coef))
    if key in seen: continue
    seen.add(key)
    if all(r==p or (F.split(r) and r<Y) for r,e in fac) and any(r!=p for r,e in fac): smooth.append((c,coef,fac))
fb=sorted(set(r for c,coef,fac in smooth for r,e in fac if r!=p)); pos={q:i for i,q in enumerate(fb)}
print("distinct elements: %d, smooth over split primes < %d: %d, primes that occur: %d: %s"%(len(seen),Y,len(smooth),len(fb),fb[:20]))
R={q:F.periods_mod(q) for q in fb}; mat=[]; used=0
for c,coef,fac in smooth:
    pats={}; ok=True
    for r,e in fac:
        if r==p: continue
        pat=[1 if (c+sum(a*R[r][(i+j)%N] for i,a in coef))%r==0 else 0 for j in range(N)]
        if sum(pat)!=e: ok=False; break
        pats[r]=pat
    if not ok: continue
    used+=1
    for s in range(N):
        v=[0]*(N*len(fb))
        for r,pat in pats.items():
            for j in range(N): v[pos[r]*N+j]=pat[(j+s)%N]
        mat.append(v)
print("elements with valuations 0/1: %d; vectors: %d; dimension: %d"%(used,len(mat),N*len(fb)))
M=flint.fmpz_mat(len(mat),N*len(fb),[x for v in mat for x in v]); rk=M.rank(); print("rank",rk)
if rk==N*len(fb):
    H=M.hnf(); d=1
    for i in range(N*len(fb)): d*=int(H[i,i])
    print("index of the lattice of relations:",d)
