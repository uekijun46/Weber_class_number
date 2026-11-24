# Second, independent check of a file of rows in the range format (relations/README.md).  It does not import verify.py
# and uses PARI instead of FLINT:
#   - norm:       |Res(Phi_M, A)| = |N(alpha)|^(p-1), where A = c + sum_i a_i sum_{h in H} x^(g^i h mod M)   (PARI polresultant)
#                 For n = 1 the resultant is slow (degree p(p-1)), and the norm is computed as the determinant of the matrix of
#                 the multiplication by alpha on the basis 1, eta_0, ..., eta_(p-2) of K (exact integers, PARI matdet); the
#                 resultant is still used for one row in every 2000.  The products eta_i * eta_j are obtained by counting
#                 exponents: sum_{a,b in H} zeta^(g^i a + g^j b), where the exponents prime to p form whole cosets of H, the sum
#                 of all N periods is 0, and the exponents divisible by p give 1 (exponent 0) or the p-th roots of unity != 1,
#                 on which H acts simply transitively (their sum is -1).
#   - primality:  PARI isprime (a proof, not a pseudoprimality test)
#   - the list of the split primes in [A,B): every q in the residue classes with q^(p-1) = 1 mod M, tested by isprime
# Checked for every row "q | c | i:a ... | r^e ...":  the coefficients are well formed; every r is prime and the r are distinct;
# the norm equals prod r^e; q is a split prime with A <= q < B and exponent 1; every other split prime factor is < q.
# Checked for the file: every split prime in [A,B) has a row.
#   python3 relations_recheck.py p n FILE A B [processes]
import sys, json, hashlib
from multiprocessing import Pool
p,n,fn,A,B=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],int(sys.argv[4]),int(sys.argv[5])
procs=int(sys.argv[6]) if len(sys.argv)>6 else 1
N=p**n; M=p**(n+1); phi=(p-1)*N
def order(x):
    k=1; y=x%M
    while y!=1: y=y*x%M; k+=1
    return k
g=next(x for x in range(2,M) if x%p and order(x)==phi)             # least generator of (Z/M)^*
H=[pow(g,N*t,M) for t in range(p-1)]
E=[[pow(g,i,M)*h%M for h in H] for i in range(N)]
if sorted(e for ex in E for e in ex)!=[x for x in range(1,M) if x%p]: raise SystemExit("FAILED: cosets")
def split(q): return q%p!=0 and pow(q,p-1,M)==1
MULT=None
if n==1:                                                            # matrices of the multiplication by eta_i on 1, eta_0, ..., eta_(N-2)
    coset={e:i for i in range(N) for e in E[i]}
    def coords(cnt):                                                # cnt: exponent mod M -> multiplicity (an H-invariant multiset)
        m=[0]*N; zero=cnt.get(0,0); pth=0
        for e,k in cnt.items():
            if e%p: m[coset[e]]+=k
            elif e: pth+=k
        for i in range(N):
            if any(cnt.get(e,0)*(p-1)!=m[i] for e in E[i]): raise SystemExit("FAILED: multiset is not invariant under H")
        if pth%(p-1): raise SystemExit("FAILED: p-th roots of unity")
        m=[x//(p-1) for x in m]
        return [zero-pth//(p-1)]+[m[k]-m[N-1] for k in range(N-1)]  # eta_(N-1) = -(eta_0 + ... + eta_(N-2))
    def period_product(i,j):
        cnt={}
        for a in E[i]:
            for b in E[j]: cnt[(a+b)%M]=cnt.get((a+b)%M,0)+1
        return coords(cnt)
    unit=lambda k:[1 if t==k else 0 for t in range(N)]
    etac=[unit(k+1) for k in range(N-1)]+[[0]+[-1]*(N-1)]          # coordinates of eta_0, ..., eta_(N-1)
    # MULT[i][k] = coordinates of eta_i * (k-th basis element)
    MULT=[[etac[i]]+[period_product(i,j) for j in range(N-1)] for i in range(N)]
_pari=None
def P():                                                           # PARI is started in the worker processes only
    global _pari,_Phi
    if _pari is None:
        import cypari2
        _pari=cypari2.Pari(); _pari.allocatemem(10**9); _Phi=_pari("polcyclo(%d)"%M)
    return _pari
def check(line):
    pari=P(); errs=[]
    try:
        qs,cs,co,fa=line.split("|"); q=int(qs); c=int(cs)
        coef=[tuple(int(u) for u in t.split(":")) for t in co.split()]
        fac=[tuple(int(u) for u in t.split("^")) for t in fa.split()]
        if any(len(t)!=2 for t in coef+fac): raise ValueError
    except ValueError:
        return None,["cannot read the row: "+line[:60]]
    idx=[i for i,a in coef]
    if len(set(idx))!=len(idx) or any(not 0<=i<N for i in idx) or any(a==0 for i,a in coef): errs.append("coefficients")
    if len(set(r for r,e in fac))!=len(fac): errs.append("repeated prime")
    prod=1
    for r,e in fac:
        if e<1 or r<2 or not bool(pari.isprime(r)): errs.append("factor %d is not prime"%r)
        prod*=r**e
    if not errs and MULT is not None:
        rows=[[(c if k==t else 0)+sum(a*MULT[i][k][t] for i,a in coef) for t in range(N)] for k in range(N)]
        if abs(int(pari.matrix(N,N,[x for r in rows for x in r]).matdet()))!=prod: errs.append("norm")
    if not errs and (MULT is None or q%2000<2):
        v=[0]*M; v[0]=c
        for i,a in coef:
            for e in E[i]: v[e]+=a
        R=int(pari.polresultant(_Phi,pari("Pol(%s)"%str(v[::-1]).replace(" ",""))))
        if abs(R)!=prod**(p-1): errs.append("norm (resultant)")
    if not (q>1 and bool(pari.isprime(q)) and split(q) and dict(fac).get(q)==1): errs.append("q is not a split prime of exponent 1")
    if not A<=q<B: errs.append("outside the range")
    for r,e in fac:
        if r!=q and r!=p and split(r) and not r<q: errs.append("split prime factor %d is not smaller than q"%r)
    return q,errs
def split_primes(job):
    pari=P(); r,lo,hi=job; q=lo-(lo%M)+r
    if q<lo: q+=M
    out=[]
    while q<hi:
        if q>1 and bool(pari.isprime(q)): out.append(q)
        q+=M
    return out
if __name__=="__main__":
    data=open(fn,"rb").read(); sha=hashlib.sha256(data).hexdigest()
    text=data.decode().split("\n")
    if text[0].split()!=["#","field",str(p),str(n)]: raise SystemExit("FAILED: the first line must be '# field %d %d'"%(p,n))
    lines=[l for l in (t.strip() for t in text) if l and not l.startswith("#")]      # as in verify.py
    step=max(1,(B-A)//(4*procs)+1); jobs=[(r,lo,min(B,lo+step)) for r in range(1,M) if split(r) for lo in range(A,B,step)]
    with Pool(procs) as pool:
        res=pool.map(check,lines,chunksize=500)
        sp=set(q for part in pool.map(split_primes,jobs) for q in part)
    bad=[(q,e) for q,e in res if e]
    for q,e in bad[:20]: print("BAD",q,e)
    qs=set(q for q,e in res if not e)
    missing=sorted(sp-qs); extra=sorted(qs-sp)
    if missing: print("split primes without a valid row: %d, e.g. %s"%(len(missing),missing[:5]))
    if extra: print("rows for numbers that are not split primes in the range: %d, e.g. %s"%(len(extra),extra[:5]))
    if bad or missing or extra: raise SystemExit("FAILED: %d bad rows, %d split primes without a row, %d extra"%(len(bad),len(missing),len(extra)))
    print("RECHECK PASS "+json.dumps({"field":[p,n],"range":[A,B],"rows":len(lines),"split_primes":len(sp),"rows_sha256":sha,
        "method":("PARI matdet of the regular representation, polresultant on a sample, isprime" if n==1 else "PARI polresultant and isprime")+"; independent of verify.py"}))
