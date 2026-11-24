# Class number 1 for the proper subfields k_{p^j} (1 <= j < n) of k_{p^n}, with no input from the literature except
# Miller's lemma, the rank theorem and Iwasawa's theorem.
#   python3 subfields.py p n
# For j = 1, 2, ..., n-1 in this order: Miller's lemma with the empty set of primes gives h(k_{p^j}) < bound.  If the bound is
# below the least value of l^(order of l mod p^j) over the primes l != p, then no prime l != p divides h(k_{p^j}) (rank theorem;
# it uses class number 1 for k_{p^i}, i < j, which was shown in the earlier steps), and p does not divide it (Iwasawa).
import sys
from verify import K, bound, up, need
p,n=int(sys.argv[1]),int(sys.argv[2]); ok=True
for j in range(1,n):
    F=K(p,j); hm=F.hmin(); best=None
    for c in [2,3,4,5,6,8,10,12,15,20,30]:
        I_up,_,B,hb=bound(F,[],c)
        if hb is not None and (best is None or hb<best[1]): best=(c,hb)
    good=best is not None and best[1]<hm; ok=ok and good
    print("k_{%d^%d}: degree %d, discriminant %d^%d, "%(p,j,F.N,p,F.logrd()[1])+("no bound" if best is None else "h < %s (c=%d)"%(up(best[1],4),best[0]))+", least possible class number > 1: %d  =>  %s"%(hm,"h = 1" if good else "NOT DECIDED"))
need(ok,"some proper subfield is not decided")
print("VERDICT: PASS  (every proper subfield of k_{%d^%d} has class number 1)"%(p,n))
