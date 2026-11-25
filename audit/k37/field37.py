#!/usr/bin/env python3
# Audit of k_37 (independent of relations/span.py and relations/verify.py; nothing is imported from relations/).
# Exact model of K = k_p (first layer): the subfield of degree p of Q(zeta), zeta of order M = p^2.
#   H = {x in (Z/M)^* : x^(p-1) = 1},  g = least primitive root modulo M,  E_i = g^i H,  eta_i = sum_{e in E_i} zeta^e  (i mod p).
#   sum_i eta_i = 0, so 1, eta_0, ..., eta_(p-2) is a Q-basis of K (checked below by a nonzero determinant).
# Products: eta_i eta_j = sum_{a in E_i, b in E_j} zeta^(a+b).  The multiset of exponents is H-invariant.
#   exponent 0: contributes 1.   exponent = p k != 0: the H-orbit is the set of all nonzero multiples of p (H -> (Z/p)^* is bijective),
#   its sum is -1.   exponent prime to p: the H-orbit is a coset E_k, its sum is eta_k.
# Everything is integer arithmetic.  Checks are explicit (no assert).
import sys
from fractions import Fraction
def need(c,msg):
    if not c: raise SystemExit("AUDIT FAILED: "+str(msg))
def is_prime_small(x):
    if x<2: return False
    d=2
    while d*d<=x:
        if x%d==0: return False
        d+=1
    return True
class Field:
    def __init__(self,p):
        need(is_prime_small(p) and p>2,"p"); self.p=p; self.M=M=p*p
        phi=p*(p-1); prs=[r for r in range(2,phi+1) if phi%r==0 and is_prime_small(r)]
        def order_full(x): return x%p!=0 and all(pow(x,phi//r,M)!=1 for r in prs)
        self.g=g=next(x for x in range(2,M) if order_full(x))
        self.H=sorted(x for x in range(1,M) if x%p and pow(x,p-1,M)==1); need(len(self.H)==p-1,"H")
        need(sorted(h%p for h in self.H)==list(range(1,p)),"H -> (Z/p)^* is not bijective")
        self.E=[sorted(pow(g,i,M)*h%M for h in self.H) for i in range(p)]
        need(sorted(e for E in self.E for e in E)==[x for x in range(M) if x%p],"cosets")
        need(sorted(pow(g,p,M)*h%M for h in self.H)==self.H,"g^p is not in H")      # so eta_(i+p) = eta_i
        self.coset={e:i for i,E in enumerate(self.E) for e in E}
        # structure constants: prod[i][j] = (constant, [m_0..m_(p-1)]) with eta_i eta_j = constant + sum_k m_k eta_k  (before using sum eta = 0)
        self.prod=[[None]*p for _ in range(p)]
        for i in range(p):
            for j in range(i,p):
                cnt=[0]*M
                for a in self.E[i]:
                    for b in self.E[j]: cnt[(a+b)%M]+=1
                const=cnt[0]
                mult_p=set(cnt[p*k] for k in range(1,p)); need(len(mult_p)==1,"multiplicity on multiples of p is not constant")
                const-=mult_p.pop()
                m=[]
                for k in range(p):
                    s=set(cnt[e] for e in self.E[k]); need(len(s)==1,"multiplicity on a coset is not constant"); m.append(s.pop())
                self.prod[i][j]=self.prod[j][i]=(const,m)
    # elements as vectors of length p in the basis 1, eta_0, ..., eta_(p-2)
    def reduce(self,const,m):                        # const + sum_{k<p} m_k eta_k  ->  basis vector (eta_(p-1) = -sum_{k<p-1} eta_k)
        p=self.p; return [const]+[m[k]-m[p-1] for k in range(p-1)]
    def elem(self,c,coef):                           # c + sum a_i eta_i
        m=[0]*self.p
        for i,a in coef: m[i%self.p]+=a
        return self.reduce(c,m)
    def mul_matrix(self,v):                          # rows: v * (basis element k), k = 0..p-1, in the basis
        p=self.p; rows=[list(v)]
        full=[v[1+k] for k in range(p-1)]+[0]        # v = v0 + sum_k full_k eta_k
        for b in range(p-1):                         # v * eta_b
            const=0; m=[0]*p; m[b]+=v[0]
            for k in range(p-1):
                if full[k]:
                    c2,m2=self.prod[k][b]; const+=full[k]*c2
                    for t in range(p): m[t]+=full[k]*m2[t]
            rows.append(self.reduce(const,m))
        return rows
    def shift(self,c,coef,j): return c,[((i+j)%self.p,a) for i,a in coef]         # tau^j
    def split(self,q): return q!=self.p and pow(q,self.p-1,self.M)==1
def read_cert(fn):
    out=[]
    for ln,line in enumerate(open(fn)):
        line=line.strip()
        if not line or line.startswith("#"): continue
        a,b,c=line.split("|")
        coef=[(int(t.split(":")[0]),int(t.split(":")[1])) for t in b.split()]
        fac=[(int(t.split("^")[0]),int(t.split("^")[1])) for t in c.split()]
        out.append((int(a),coef,fac))
    return out
if __name__=="__main__":
    F=Field(int(sys.argv[1]) if len(sys.argv)>1 else 37)
    print("p=%d M=%d g=%d |H|=%d"%(F.p,F.M,F.g,len(F.H)))
