# Write a certificate file from the rows found by search.py (untrusted helper; everything is re-checked by verify.py).
#   python3 export.py p n X seedprime
import sys, pickle, os
import flint
from field import Field
from search import peel
p,n,X,q0=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
F=Field(p,n); N=F.N
rows=pickle.load(open("rows_%d_%d.pkl"%(p,n),"rb"))
extra=pickle.load(open("targeted_%d_%d.pkl"%(p,n),"rb")) if os.path.exists("targeted_%d_%d.pkl"%(p,n)) else []
rows=rows+[(c,coef,fac) for (q,c,coef,fac) in extra]
# --- lattice block for the seed prime q0: elements with norm p^a q0^k lying in exactly k distinct primes above q0
R=F.periods_mod(q0)
pure=sorted(set((c,tuple(coef),tuple(fac)) for (c,coef,fac) in rows if all(r in (p,q0) for r,e in fac) and any(r==q0 for r,e in fac)),key=lambda t:len(t[1]))
block=[]; mat=[]
def hnf_is_identity(mat):
    Hn=flint.fmpz_mat(len(mat),N,[x for row in mat for x in row]).hnf()
    return all(int(Hn[i,j])==(1 if i==j else 0) for i in range(N) for j in range(N))
for (c,coef,fac) in pure:
    pat=[1 if (c+sum(a*R[(i+j)%N] for i,a in coef))%q0==0 else 0 for j in range(N)]
    if sum(pat)!=dict(fac)[q0]: continue                       # some valuation exceeds 1: skip
    block.append((c,coef)); mat+= [pat[s:]+pat[:s] for s in range(N)]
    if len(block)>=3 and hnf_is_identity(mat): break
assert hnf_is_identity(mat), "lattice block for the seed prime not found"
print("lattice block for %d: %d elements"%(q0,len(block)))
# --- rows needed for the split primes <= X, with their dependencies, in an order in which every dependency comes first
cert,order=peel(F,rows,[q0])
sv=bytearray([1])*(X+1)
for i in range(2,int(X**.5)+1):
    if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
targets=[q for q in range(2,X+1) if sv[q] and F.split(q)]
assert all(q in cert for q in targets), "some target primes are not certified"
needed=set(); stack=[q for q in targets if q!=q0]
while stack:
    q=stack.pop()
    if q in needed or q==q0: continue
    needed.add(q)
    for r,e in rows[cert[q]][2]:
        if r!=q and F.split(r) and r not in needed: stack.append(r)
out=open("cert_k%d_%d.txt"%(p,n),"w")
out.write("# field %d %d\n"%(p,n))
out.write("L %d | %s\n"%(q0," ; ".join("%d,%s"%(c," ".join("%d:%d"%ia for ia in coef)) for c,coef in block)))
cnt=0
for q in order:
    if q in needed:
        c,coef,fac=rows[cert[q]]
        out.write("%d | %d | %s | %s\n"%(q,c," ".join("%d:%d"%ia for ia in coef)," ".join("%d^%d"%re for re in fac))); cnt+=1
out.close()
print("rows written: %d (targets %d, helper primes above X: %d); file size %.1f MB"%(cnt,len(targets),sum(1 for q in needed if q>X),os.path.getsize("cert_k%d_%d.txt"%(p,n))/1e6))
