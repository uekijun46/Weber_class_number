#!/usr/bin/env python3
# Audit k_47: the 2 x 2 block of the two primes b = 6453089, B = 184398551, which occur exactly in the data lines 5417 and 10574.
import json
p=47; raw=json.load(open("out/rows_47.json")); b,B=6453089,184398551; L1,L2=5417,10574
occ={q:[i for i,r in enumerate(raw) if any(q2==q for q2,J in r)] for q in (b,B)}
print("lines containing b: %s; lines containing B: %s"%(occ[b],occ[B]))
e={(i,q):dict((q2,J) for q2,J in raw[i])[q] for i in (L1,L2) for q in (b,B)}
print("entries (sets J; the entry is sum_{j in J} x^j):",{k:v for k,v in e.items()})
a1,c1,a2,c2=e[(L1,b)][0],e[(L1,B)][0],e[(L2,b)][0],e[(L2,B)][0]
m=(a1+c2-a2-c1)%p
print("block [[x^%d, x^%d],[x^%d, x^%d]]: det = x^%d - x^%d = x^%d (1 - x^%d) up to sign; m = %d mod 47 (m != 0: %s)"%(a1,c1,a2,c2,(a1+c2)%p,(a2+c1)%p,(a2+c1)%p,m,m,m!=0))
# Nm(1 - x^m) = Phi_47(1) = 47 for m != 0 mod 47; check with the integer matrix of multiplication by 1 - x^m on A
from flint import fmpz_mat
def mulx(v): top=v[-1]; return [-top]+[v[k-1]-top for k in range(1,p-1)]
rows=[]; cur=[1]+[0]*(p-2)
for j in range(p-1):
    y=cur[:]
    for _ in range(m): y=mulx(y)
    rows.append([u-w for u,w in zip(cur,y)]); cur=mulx(cur)
print("Nm(1 - x^m) by the determinant of the multiplication matrix:",int(fmpz_mat(p-1,p-1,[c for r in rows for c in r]).det()))
print("other primes of the two lines:",[q for q,J in raw[L1] if q not in (b,B)],[q for q,J in raw[L2] if q not in (b,B)])
