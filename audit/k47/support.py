#!/usr/bin/env python3
# Audit k_47: the supports of the columns (which lines of the certificate contain which prime), from the certificate alone.
# Looks for pairs of primes that occur in exactly the same two lines (the author's guess about the origin of the factor 47).
#   python3 support.py
import collections
from fieldp import read_cert, CERT
p=47; cert=read_cert(CERT(p),p)
occ=collections.defaultdict(list)
for n,(c,coef,fac) in enumerate(cert):
    for q,e in fac:
        if q!=p: occ[q].append(n)
cnt=collections.Counter(len(v) for v in occ.values())
print("lines %d; primes %d; number of primes by the number of lines in which they occur: %s ..."%(len(cert),len(occ),sorted(cnt.items())[:8]))
big=collections.Counter(sum(1 for q,e in fac if q>=10**6) for c,coef,fac in cert)
print("number of prime factors >= 10^6 per line:",sorted(big.items()))
bysup=collections.defaultdict(list)
for q,v in occ.items(): bysup[tuple(v)].append(q)
same=[(sup,sorted(qs)) for sup,qs in bysup.items() if len(qs)>1]
print("groups of primes with the same support:",len(same))
for sup,qs in sorted(same,key=lambda t:t[1]):
    print("  primes %s occur exactly in the lines %s (0-based data lines)"%(qs,list(sup)))
    for n in sup: print("     line %d: %s | %s | %s"%(n,cert[n][0]," ".join("%d:%d"%t for t in cert[n][1])," ".join("%d^%d"%t for t in cert[n][2])))
