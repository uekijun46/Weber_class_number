#!/usr/bin/env python3
# Audit k_47: tests of the changed acceptance rule of relations/span.py (imports span.py on purpose: it is the object under test).
#  (1) span.ppower against a direct definition.
#  (2) span.prove on small synthetic rows over A = Z[x]/(Phi_47) with known index; the acceptance test of span.main
#      (g != 0 and ppower(g, p)) is applied to the returned gcd.   (prove only uses F.N.)
#   python3 ppower_test.py
import sys, os, types, io, contextlib
here=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(here,"..","..","relations"))
import span
from flint import fmpz_poly
def need(c,msg):
    if not c: raise SystemExit("AUDIT FAILED: "+str(msg))
p=47
pw={p**e for e in range(0,40)}
for g in list(range(1,20000))+[p**k*m for k in range(0,6) for m in (1,2,3,46,48,283,2209+1)]+[p**30,p**30*2,p**30+1]:
    need(span.ppower(g,p)==(g in pw),"ppower(%d)"%g)
print("(1) span.ppower(g, 47) agrees with the definition for g = 1..19999 and 45 further values (ppower(0, p) would not terminate; span.py tests g != 0 first in all three places)")
F=types.SimpleNamespace(N=p)
one=fmpz_poly([1]); omx=fmpz_poly([1,-1]); x=fmpz_poly([0,1])
def run(name,rows,core,expect_g,expect_accept):
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf): T,g=span.prove(F,[dict(r) for r in rows],core)
    acc=(g!=0 and span.ppower(g,p))
    print("  %-58s gcd = %s; accepted by span.main's rule: %s"%(name,g if g<10**12 else "%d digits"%len(str(g)),acc))
    need((expect_g is None or g==expect_g) and acc==expect_accept,name)
print("(2) span.prove on synthetic rows (1 or 2 columns q = 101, 103; core = number of columns, so that step (c) decides):")
# a) identity-like: generates everything
run("a) rows e1, e2 (twice each)",[{101:one},{103:one},{101:one},{103:one}],2,1,True)
# b) the situation of k_47: two columns that occur in the same two rows with unit entries
run("b) [x^23, x^9], [x^13, x^28] (twice)",[{101:x**23,103:x**9},{101:x**13,103:x**28}]*2,2,47,True)
# c) index 47^2
run("c) [1-x, 0], [0, 1-x] (twice)",[{101:omx},{103:omx}]*2,2,47**2,True)
# d) index Nm(2)^2 = 2^92: must be rejected   (the sets of step (c) are consecutive rows in the order of the trace form, so rows of equal size are used)
run("d) [2, 0], [0, 2] (twice)",[{101:fmpz_poly([2])},{103:fmpz_poly([2])}]*2,2,2**92,False)
# e) index (47 * 2^46)^2: must be rejected
run("e) [2(1-x), 0], [0, 2(1-x)] (twice)",[{101:2*omx},{103:2*omx}]*2,2,(47*2**46)**2,False)
# f) rank deficient: every determinant is 0, gcd 0: must be rejected (and must not hang)
run("f) [1-x, 1-x], [1-x, 1-x] (twice)",[{101:omx,103:omx},{101:omx,103:omx}]*2,2,0,False)
# g) Nm(x - 2) = 2^47 - 1 = 2351 * 4513 * 13264529 (prime factors = 1 mod 47): must be rejected
run("g) [x-2, 0], [0, x-2] (twice)",[{101:fmpz_poly([-2,1])},{103:fmpz_poly([-2,1])}]*2,2,(2**47-1)**2,False)
# h) gcd of two different multiples: [2(1-x)] and [3(1-x)] on one column -> gcd(47*2^46, 47*3^46) = 47: accepted, correctly (the two rows together generate (1-x))
run("h) one column, rows 2(1-x) and 3(1-x)",[{101:2*omx},{101:3*omx}],1,47,True)
print("ppower_test passed")
