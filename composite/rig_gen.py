# Certificate GENERATOR (untrusted; uses PARI bnfinit/bnfisprincipal as an oracle).  Output is checked by rig_verify.py.
#
# Usage:  python3 rig_gen.py FIELD_INDEX X [TIME_LIMIT_SECONDS]
#   FIELD_INDEX = index into rig_common.FIELDS (0..12) or a name like 2_1x13_1.
#   Writes / extends rig_cert_<name>.txt :
#     line 1 : T  (integer coefficients, low degree first)            -- polredbest of the Gaussian-period polynomial
#     line 2 : #EMBED D ; b_0 ... b_{N-1}   meaning  theta = (sum b_k eta^k)/D  is a root of T, eta = Gaussian period
#     then   : q ; D ; a_0 ... a_{N-1}      meaning  alpha = (sum a_k theta^k)/D has norm +-q (q totally split, or a
#              ramified prime of residue degree 1 that splits completely in the other layers)
#   Resumable: rerun with the same or larger X; primes already present are skipped.
#   Also writes rig_units_<name>.txt (units of the subfield K' for the Chevalley check) when it is needed.
import sys, os, time, math
import cypari2
from rig_common import *

pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9); pari.default("realprecision", 100)

def polstr(P): return " ".join(str(int(c)) for c in P.coeffs())
def topari(P): return pari("Pol(%s)" % ("[" + ",".join(str(int(c)) for c in reversed(P.coeffs())) + "]"))
def coeffs_with_den(a, N):
    """a = PARI polynomial (t_POL, rational coeffs) of degree < N -> (D, [ints])"""
    a = pari.lift(a); D = int(pari.denominator(pari.content(a))); v = pari.Vecrev(a * D, N)
    return D, [int(c) for c in v]

def get_field(F):
    nm = name(F); N = degree(F); fn = "rig_cert_%s.txt" % nm
    if os.path.exists(fn):
        lines = open(fn).read().split("\n")
        T = pari("Pol([%s])" % ",".join(reversed(lines[0].split())))
        return fn, T, None
    T0, Z, eta, _ = period_minpoly(F)
    R = pari.polredbest(topari(T0), 1); T = R[0]
    b = pari.modreverse(R[1])                 # Mod(b(x), T0): a root of T as a polynomial in a root (eta) of T0
    D, bc = coeffs_with_den(b, N)
    with open(fn, "w") as f:
        f.write(" ".join(str(int(c)) for c in pari.Vecrev(T)) + "\n")
        f.write("#EMBED %d ; %s\n" % (D, " ".join(map(str, bc))))
    return fn, T, None

def done_primes(fn):
    s = set()
    for line in open(fn):
        if line.startswith("#") or ";" not in line: continue
        s.add(int(line.split(";")[0]))
    return s

def need_cert(q, F):
    e, f, g = decomposition(q, F)
    return f == 1 and not free_by_subfield(e, f)

def units_file(F):
    """For each ramified p_i with >1 prime in K' (compositum of the other layers): units of K' in the eta'-basis."""
    out = []
    for p, n in F:
        Fo = [(q, k) for q, k in F if q != p]
        e, f, g = decomposition(p, F)                # g = number of primes of K above p = number in K'
        if g == 1: continue
        T0, Z, eta, _ = period_minpoly(Fo); Np = degree(Fo)
        bnf = pari.bnfinit(topari(T0), 1)
        us = list(bnf.bnf_get_fu())
        out.append("#KPRIME %d %s\n" % (p, repr(Fo)))
        for u in us:
            D, c = coeffs_with_den(pari.lift(u), Np)
            out.append("%d ; %s\n" % (D, " ".join(map(str, c))))
    if out:
        open("rig_units_%s.txt" % name(F), "w").write("".join(out))

def main():
    arg = sys.argv[1]
    F = FIELDS[int(arg)] if arg.isdigit() else [tuple(map(int, s.split("_"))) for s in arg.split("x")]
    X = int(float(sys.argv[2])); LIM = float(sys.argv[3]) if len(sys.argv) > 3 else 500
    t0 = time.time(); N = degree(F)
    fn, T, _ = get_field(F)
    units_file(F)
    have = done_primes(fn)
    todo = [q for q in small_primes(X) if q not in have and need_cert(q, F)]
    print(name(F), "N=%d" % N, "to do:", len(todo), "have:", len(have), flush=True)
    if not todo: return
    bnf = pari.bnfinit(T, 1)
    print("bnfinit %.1fs" % (time.time() - t0), flush=True)
    out = open(fn, "a"); n = 0
    for q in todo:
        if time.time() - t0 > LIM: print("time limit reached at q =", q); break
        P = pari.idealprimedec(bnf, q)
        e, f, g = decomposition(q, F); assert len(P) == g
        r = pari.bnfisprincipal(bnf, P[0])
        if any(int(x) for x in r[0]): print("NON-PRINCIPAL prime above", q); break
        a = pari.nfbasistoalg(bnf, r[1])
        assert abs(int(pari.norm(a))) == q
        D, c = coeffs_with_den(a, N)
        out.write("%d ; %d ; %s\n" % (q, D, " ".join(map(str, c)))); out.flush(); n += 1
    out.close()
    print("wrote %d certificates, %.0fs" % (n, time.time() - t0))

if __name__ == "__main__": main()
