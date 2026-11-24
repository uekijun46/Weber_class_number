# Bootstrapping the smallest split prime q = 19 of k_{7^2} without any class group computation:
# for elements alpha whose norm is 7^a * 19^k, compute the 19-adic valuations of alpha at the 49 primes above 19;
# the cyclic shifts of these vectors are relations among the classes of the 49 primes.  If they span Z^49, the primes are principal.
import pickle, cypari2, time
from proto import split, NN, M
pari=cypari2.Pari(); pari.allocatemem(10**9,silent=True)
from mpmath import mp, mpf, cos, pi
mp.prec=400
def tomp(x): return mp.mpf(str(x).replace(' ',''))
pari.default('realprecision',100)
g=next(x for x in range(2,M) if x%7 and pow(x,147,M)!=1 and pow(x,98,M)!=1 and pow(x,42,M)!=1)
H=[pow(g,NN*t,M) for t in range(6)]
eta=[sum(cos(2*pi*(pow(g,i,M)*h%M)/M) for h in H) for i in range(NN)]
q=19; k=8; Q=q**k
# the 49 embeddings O_K -> Z_19: eta_i -> sum_{h in H} z^(g^(i+j) h), z a primitive 343rd root of unity in the unramified
# extension of Z_19 of degree 6 (the order of 19 modulo 343), computed modulo 19^k
fac=pari.factorpadic(pari.polcyclo(343),q,k)[0]; gpol=pari.lift(fac[0]); assert int(pari.poldegree(gpol))==6
z=pari.Mod(pari("y")*pari.Mod(1,Q), pari.subst(gpol,"x","y")*pari.Mod(1,Q))
R=[]
for j in range(NN):
    v=sum(z**(pow(g,j,M)*h%M) for h in H); v=pari.lift(pari.lift(v))
    assert int(pari.poldegree(v))<=0, "period not in Z_19"
    R.append(int(pari.polcoef(v,0))%Q)
assert all(any(R[(i+t)%NN]%q!=R[i]%q for i in range(NN)) for t in range(1,NN))   # the 49 shifted sequences are distinct mod 19; print(len(set(x%q for x in R)))
rows=pickle.load(open("rows.pkl","rb"))
pure=sorted(set((c,idx,sg,fac) for (c,idx,sg,fac) in rows if set(p for p,e in fac if split(p))=={19} and all(p in (7,19) for p,e in fac)))
print("distinct elements with norm 7^a * 19^k:",len(pure))
def val(x):
    v=0
    while x%q==0 and v<k: x//=q; v+=1
    return v
vecs=[]
for (c,idx,sg,fac) in pure:
    v=[val((c+sum(s*R[(i+j)%NN] for i,s in zip(idx,sg)))%Q) for j in range(NN)]
    assert sum(v)==dict(fac)[19] and max(v)<k, (v,fac)
    vecs.append(v)
import collections
print("valuation patterns (max valuation at one prime -> count):",sorted(collections.Counter(max(v) for v in vecs).items()))
def rank_lattice(vs):
    mat=[]
    for v in vs:
        for t in range(NN): mat.append(v[t:]+v[:t])
    Mx=pari.matrix(len(mat),NN,[x for row in mat for x in row])
    Hn=pari.mathnf(pari.mattranspose(Mx))
    return int(pari.matrank(Mx)), abs(int(pari.matdet(Hn))) if int(pari.matsize(Hn)[1])==NN else None
for n in (1,2,3,5,10,len(vecs)):
    rk,idx=rank_lattice(vecs[:n]); print("using %d elements: rank %d, index of the relation lattice in Z^49: %s"%(n,rk,idx))
    if idx==1: break
print("first elements used:",[(c,idx,sg,fac) for (c,idx,sg,fac) in pure[:3]])
