# Independent RIGOROUS verifier.  No PARI: only exact integer / polynomial arithmetic (Python ints, python-flint
# fmpz/fmpz_poly/fmpz_mat/nmod_poly) and ball arithmetic (arb/acb, with acb.integral for the integral).
#
# Usage:  python3 rig_verify.py NAME_OR_INDEX [--X XMAX] [--cs 6,8,10,...] [--qfree QMAX]
#   NAME like 2_1x13_1 or an index into rig_common.FIELDS; reads rig_cert_<name>.txt (and rig_units_<name>.txt).
#   --X     only use certificates for primes <= XMAX (to compare with an earlier run); default: all.
#   Prints the checks performed and a final line  RESULT <name> ... PASS/FAIL.
import sys, math, time
from flint import fmpz, fmpz_poly, fmpz_mat, nmod_poly, arb, acb, ctx
from rig_common import *

def log(*a): print(*a, flush=True)

# ---------------------------------------------------------------- field
def read_cert(F):
    fn = "rig_cert_%s.txt" % name(F); lines = [l for l in open(fn).read().split("\n") if l.strip()]
    T = fmpz_poly([int(x) for x in lines[0].split()])
    assert lines[1].startswith("#EMBED")
    Ds, bs = lines[1][6:].split(";"); emb = (int(Ds), [int(x) for x in bs.split()])
    certs = []
    for l in lines[2:]:
        if l.startswith("#"): continue
        q, D, a = l.split(";"); certs.append((int(q), int(D), [int(x) for x in a.split()]))
    return T, emb, certs

def check_field(F, T, emb):
    N = degree(F); m = conductor(F)
    T0, Z, eta, _ = period_minpoly(F)                       # exact: T0(eta)=0, conjugates distinct => T0 = minpoly(eta)
    log("  K = fixed field of H in Q(zeta_%d), |H| = %d; eta = Gaussian period, minpoly T0 of degree %d verified" %
        (m, len(group_H(F)), T0.degree()))
    assert T.degree() == N and T.leading_coefficient() == 1
    c, fac = T.factor(); assert len(fac) == 1 and fac[0][1] == 1, "T not irreducible"
    D, b = emb; assert len(b) <= N and D > 0
    y = eval_rational_at(Z, b, D, eta)                      # y = D*theta in Z[zeta_m]
    assert horner_exact(Z, T, y, D) == 0, "theta is not a root of T"
    log("  T irreducible of degree %d, monic, and T(theta)=0 exactly for theta = b(eta)/%d in Q(zeta_%d): Q(theta) = K" % (N, D, m))
    # cross-check of the congruence data against the factorisation pattern of T mod q (not needed for the proof)
    dT = int(abs(T.discriminant())); bad = 0; cnt = 0
    for q in small_primes(3000):
        if dT % q == 0: continue
        e, f, g = decomposition(q, F); assert e == 1
        degs = []
        for P, k in nmod_poly([int(x) for x in T.coeffs()], q).factor()[1]: degs += [P.degree()] * k
        cnt += 1
        if sorted(degs) != [f] * g: bad += 1
    assert bad == 0
    log("  cross-check: factorisation pattern of T mod q equals (g,f) from congruences for all %d primes q<3000, q !| disc T" % cnt)

# ---------------------------------------------------------------- certificates
def mult_matrix(T, A):
    N = T.degree(); rows = []; v = A % T; x = fmpz_poly([0, 1])
    for i in range(N):
        c = [int(t) for t in v.coeffs()]; rows.append(c + [0] * (N - len(c)))
        v = (v * x) % T
    return fmpz_mat(rows)

def check_cert(F, T, q, D, a):
    N = T.degree()
    assert fmpz(q).is_prime(), q
    e, f, g = decomposition(q, F)
    assert f == 1 and not free_by_subfield(e, f), "certificate for a prime that is not of residue degree 1"
    assert len(a) <= N and D > 0
    M = mult_matrix(T, fmpz_poly(a))                       # multiplication by D*alpha
    cp = [int(t) for t in M.charpoly().coeffs()]            # char poly of D*alpha, exact
    assert len(cp) == N + 1 and cp[N] == 1
    Dk = 1
    for k in range(N, -1, -1):                              # charpoly of alpha: coefficient c_k / D^{N-k} must be integral
        assert cp[k] % Dk == 0, "alpha not integral"; Dk *= D
    nrm = cp[0] // D ** N * (-1) ** N
    assert abs(nrm) == q ** f, "norm check failed"
    return e, f, g

# ---------------------------------------------------------------- l | N part (Chevalley)
def hensel_factors(P, p, k):
    """monic integer lifts mod p^k of the irreducible factors of P mod p (P monic, squarefree mod p)."""
    Pm = nmod_poly([int(x) for x in P.coeffs()], p)
    assert Pm.gcd(Pm.derivative()).degree() == 0, "T' not squarefree mod p"
    facs = [h for h, mult in Pm.factor()[1]]
    out = []
    for h in facs:
        g = Pm // h; assert (g * h) == Pm
        G, s, r = h.xgcd(g); assert G.degree() == 0
        inv = pow(int(G.coeffs()[0]), -1, p); s *= inv; r *= inv          # s h + r g = 1 mod p
        H = fmpz_poly([int(x) for x in h.coeffs()]); Gz = fmpz_poly([int(x) for x in g.coeffs()]); pt = p
        for t in range(1, k):
            E = P - H * Gz
            assert all(int(c) % pt == 0 for c in E.coeffs())
            e = nmod_poly([int(c) // pt % p for c in E.coeffs()], p)
            er = e * r; qq, dh = divmod(er, h); dg = e * s + qq * g
            H = H + pt * fmpz_poly([int(c) for c in dh.coeffs()])
            Gz = Gz + pt * fmpz_poly([int(c) for c in dg.coeffs()])
            pt *= p
        E = P - H * Gz; assert all(int(c) % pt == 0 for c in E.coeffs()); assert H.leading_coefficient() == 1
        out.append(H)
    return out

def unit_check(Tp, D, a):
    N = Tp.degree(); M = mult_matrix(Tp, fmpz_poly(a)); cp = [int(t) for t in M.charpoly().coeffs()]
    Dk = 1
    for k in range(N, -1, -1): assert cp[k] % Dk == 0; Dk *= D
    assert abs(cp[0]) == D ** N, "not a unit"

def chevalley(F):
    """For each ramified p_i: if p_i has one prime in K' (= compositum of the other layers) -> p_i !| h(K) since
       K/K' is a cyclic p-extension with a single (totally) ramified prime and p_i !| h(K').  Otherwise Chevalley:
       |Cl(K)^G| = p^{n(g'-1)} / [E_K' : E_K' cap N K^*], and the index is >= the order of the image of the given units
       under the local norm-residue maps at the g' primes above p (Hasse norm theorem: K/K' cyclic)."""
    res = {}
    for p, n in F:
        e, f, g = decomposition(p, F)
        Fo = [(q, k) for q, k in F if q != p]
        if g == 1:
            res[p] = "one prime above %d in K' (f=%d): %d !| h" % (p, f, p); continue
        lines = [l for l in open("rig_units_%s.txt" % name(F)).read().split("\n") if l.strip()]
        i = lines.index(next(l for l in lines if l.startswith("#KPRIME %d " % p)))
        us = []
        for l in lines[i + 1:]:
            if l.startswith("#"): break
            D, a = l.split(";"); us.append((int(D), [int(x) for x in a.split()]))
        Tp, Zp, etap, _ = period_minpoly(Fo)                    # K' = Q(eta'), exact
        k = n + 1 if p != 2 else n + 2; pk = p ** k; pn = p ** n
        facs = hensel_factors(Tp, p, k)
        assert len(facs) == g and all(h.degree() == f for h in facs)   # places of K' over p <-> p-adic factors
        rows = []
        for D, a in us:
            unit_check(Tp, D, a); assert D % p
            A = fmpz_poly(a); row = []
            for h in facs:
                Nv = int(h.resultant(A)) * pow(D, -h.degree(), pk) % pk   # local norm N_{K'_v/Q_p}(u) mod p^k
                assert Nv % p
                if p != 2:
                    x = pow(Nv, p - 1, pk); gen = 1 + p                    # (1+pZ)/(1+p^{n+1}Z) = <1+p>, order p^n
                else:
                    x = Nv if Nv % 4 == 1 else (-Nv) % pk; gen = 5          # Z_2^*/(+-1)(1+2^{n+2}) = <5>, order 2^n
                j = next(j for j in range(pn) if pow(gen, j, pk) == x)
                row.append(j)
            rows.append(row)
        L = fmpz_mat(rows + [[pn if i2 == j2 else 0 for j2 in range(g)] for i2 in range(g)]).hnf()
        Lr = [r for r in L.tolist() if any(r)]
        assert len(Lr) == g
        det = abs(int(fmpz_mat(Lr).det())); img = pn ** g // det
        ok = img == pn ** (g - 1)
        res[p] = "%d primes above %d in K'=%s; %d units verified; local symbols %s; image order %d (need %d): %s" % (
            g, p, name(Fo), len(us), rows, img, pn ** (g - 1), "p-part trivial" if ok else "FAILED")
        assert ok
    return res

# ---------------------------------------------------------------- analytic part
def integral_upper(c):
    c = arb(c)
    def f(x, analytic):
        Fx = (-(x / c) ** 2).exp() / (x / 2).cosh()
        return (1 - Fx) / 2 * (1 / (x / 2).sinh() + 1 / (x / 2).cosh())
    eps = arb(1) / 1000; Mx = arb(120)
    mid = acb.integral(f, acb(eps), acb(Mx)).real
    assert mid.rad() < 1e-20, "integral not accurate"
    # 0<x<eps:  1-F(x) <= x^2 (1/c^2+1/8), 1/sinh(x/2) <= 2/x, 1/cosh <= 1
    head = (1 / c ** 2 + arb(1) / 8) * (eps ** 2 / 2 + eps ** 3 / 6)
    # x>Mx: (1-F)/2 <= 1/2, 1/sinh(x/2) <= 2 e^{-x/2}/(1-e^{-Mx}), 1/cosh(x/2) <= 2 e^{-x/2}
    tail = (1 / (1 - (-Mx).exp()) + 1) * 2 * (-Mx / 2).exp()
    return arb((mid + head + tail).upper())

def prime_sum_lower(c, N, plist):
    """(2/N) sum over (q, f, g) in plist of g * sum_m f log q * q^{-fm/2} F(m f log q), truncated (all terms > 0)."""
    c = arb(c); S = arb(0)
    for q, f, g in plist:
        lq = f * arb(q).log(); qh = arb(q) ** (arb(f) / 2); s = arb(0); pw = qh
        cut = min(10 * float(c), 120.0)
        for m in range(1, 400):
            x = m * lq
            if float(x.mid()) > cut: break
            s += lq / pw * (-(x / c) ** 2).exp() / (x / 2).cosh()
            pw *= qh
        S += g * s
    return arb((2 * S / N).lower())

def main():
    args = sys.argv[1:]; arg = args[0]
    F = FIELDS[int(arg)] if arg.isdigit() else [tuple(map(int, s.split("_"))) for s in arg.split("x")]
    Xmax = int(float(args[args.index("--X") + 1])) if "--X" in args else None
    cs = [float(x) for x in args[args.index("--cs") + 1].split(",")] if "--cs" in args else [4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100]
    QF = int(float(args[args.index("--qfree") + 1])) if "--qfree" in args else 300000
    ctx.prec = 128; t0 = time.time()
    N = degree(F); m = conductor(F); nm = name(F)
    log("== %s : N=%d, conductor %d" % (nm, N, m))
    # induction order: every proper composite subfield must occur earlier in FIELDS
    idx = FIELDS.index(F) if F in FIELDS else None
    import itertools
    subs = [[(p, k) for (p, n), k in zip(F, ks) if k > 0] for ks in itertools.product(*[range(n + 1) for p, n in F])]
    subs = [S for S in subs if len(S) >= 2 and S != F]
    for S in subs: assert S in FIELDS and FIELDS.index(S) < idx, ("subfield not treated before", S)
    log("  proper composite subfields (all treated earlier in the list):", [name(S) for S in subs] or "none")
    log("  prime-power layers assumed to have h=1:", [name([(p, k)]) for p, n in F for k in range(1, n + 1)])
    T, emb, certs = read_cert(F)
    check_field(F, T, emb)
    plist = []; ncert = 0; Xc = 0; certified = set()
    for q, D, a in certs:
        if Xmax is not None and q > Xmax: continue
        e, f, g = check_cert(F, T, q, D, a); ncert += 1; Xc = max(Xc, q); certified.add(q)
        plist.append((q, f, g))
    log("  %d certificates verified (integral, norm +-q, residue degree 1), largest q = %d  (%.0fs)" % (ncert, Xc, time.time() - t0))
    # completeness information
    miss = [q for q in small_primes(Xc) if q not in certified and decomposition(q, F)[1] == 1
            and not free_by_subfield(*decomposition(q, F)[:2])]
    log("  split / non-free primes <= %d without certificate: %d %s" % (Xc, len(miss), miss[:10]))
    # free primes (a), (b)
    nfree = 0; ram = []
    for q in small_primes(QF):
        e, f, g = decomposition(q, F)
        if q in certified: continue
        if free_by_subfield(e, f):
            plist.append((q, f, g)); nfree += 1
            if e > 1: ram.append("%d:(e,f,g)=(%d,%d,%d) free" % (q, e, f, g))
        elif e > 1: ram.append("%d:(e,f,g)=(%d,%d,%d) NOT free" % (q, e, f, g))
    for p, n in F:
        e, f, g = decomposition(p, F)
        if p in certified: ram.append("%d:(e,f,g)=(%d,%d,%d) certified" % (p, e, f, g))
    log("  free primes q <= %d used: %d ; ramified: %s" % (QF, nfree, "; ".join(ram)))
    chev = chevalley(F)
    for p in chev: log("  l=%d part: %s" % (p, chev[p]))
    hm = hmin(F)
    ex = disc_exponents(F)
    logrd = sum(arb(ex[p]) / N * arb(p).log() for p in ex)
    C0 = arb.pi() / 2 + arb.const_euler() + (8 * arb.pi()).log()
    best = None
    for c in cs:
        I = integral_upper(c); S = prime_sum_lower(c, N, plist)
        B = arb((C0 - logrd - I + S).lower())
        hb = arb((2 * arb(c) * arb.pi().sqrt() / (N * B)).upper()) if B > 0 else None
        log("  c=%-5g  I(c) <= %.10f   prime sum >= %.10f   B >= %.8f   %s" % (c, float(I), float(S), float(B),
            "h < %.4f" % float(hb) if hb is not None else "B<=0"))
        if hb is not None and (best is None or hb < best[1]): best = (c, hb)
    log("  log rd = %.8f (rd = %.3f); hmin = %d" % (float(logrd), math.exp(float(logrd)), hm))
    ok = best is not None and best[1] < hm
    log("RESULT %s N=%d cond=%d certs=%d X=%d best_c=%g bound=%.4f hmin=%d %s  (%.0fs)" % (nm, N, m, ncert, Xc,
        best[0] if best else 0, float(best[1]) if best else float("inf"), hm, "PASS" if ok else "FAIL", time.time() - t0))

if __name__ == "__main__": main()
