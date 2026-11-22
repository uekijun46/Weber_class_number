# Shared exact routines for rig_gen.py / rig_verify.py.  NO PARI here: only Python integers and python-flint.
#
# K = compositum of the layers B_{p_i,n_i} of the cyclotomic Z_p-extensions, realised as the fixed field of an
# explicit subgroup H of (Z/m)^*, m = conductor.  Gal(K/Q) = (Z/m)^*/H = prod_i Gal(B_{p_i,n_i}/Q) (cyclic, orders coprime).
import math
from fractions import Fraction
from flint import fmpz, fmpz_poly, fmpq_poly, arb, arb_poly, ctx

def name(F): return "x".join("%d_%d" % (p, n) for p, n in F)

def layer_modulus(p, n): return 2 ** (n + 2) if p == 2 else p ** (n + 1)

def conductor(F): return math.prod(layer_modulus(p, n) for p, n in F)

def degree(F): return math.prod(p ** n for p, n in F)

def in_Hi(a, p, n):
    """a (coprime to p) lies in the subgroup H_i of (Z/p^{n+1})^* (resp. (Z/2^{n+2})^*) fixing B_{p,n}."""
    M = layer_modulus(p, n); a %= M
    if p == 2: return a in (1, M - 1)              # B_{2,n} = Q(zeta_{2^{n+2}})^+ : fixed by +-1
    return pow(a, p - 1, M) == 1                    # B_{p,n} = fixed field of mu_{p-1} (the unique subgroup of order p-1)

def in_H(a, F): return all(in_Hi(a, p, n) for p, n in F)

def layer_order(q, p, n):
    """order of the Frobenius of q (q != p) in Gal(B_{p,n}/Q) = (Z/M)^*/H_i  (a power of p)."""
    M = layer_modulus(p, n); k = 1
    if p == 2:
        x = q % M
        while x not in (1, M - 1): x = x * x % M; k *= 2
        return k
    x = pow(q, p - 1, M)
    while x != 1: x = pow(x, p, M); k *= p
    return k

def frob_order(q, F, skip=None):
    """residue degree f of an unramified q in the compositum of the layers of F (excluding prime `skip`).
       Gal = product of the layer groups (coprime orders), so the order of Frobenius is the product."""
    return math.prod(layer_order(q, p, n) for p, n in F if p != skip)

def decomposition(q, F):
    """(e, f, g) of the rational prime q in K, from congruences only."""
    N = degree(F); ps = [p for p, n in F]
    if q in ps:
        n = dict(F)[q]; e = q ** n                  # B_{q,n} totally ramified at q, others unramified at q
        f = frob_order(q, F, skip=q)                # residue degree in the compositum K' of the other layers
        return e, f, N // (e * f)
    f = frob_order(q, F); return 1, f, N // f

def free_by_subfield(e, f):
    """True iff every prime P (with these e,f in the cyclic K) equals P_L O_K for a prime P_L of a PROPER subfield L:
       need 1 != H' <= D (cyclic, |D| = e f) with H' cap I = 1 (|I| = e), i.e. a prime r | f with r not dividing e."""
    r = f
    for d in range(2, f + 1):
        if r % d == 0:
            if e % d: return True
            while r % d == 0: r //= d
    return False

def disc_exponents(F):
    """{p_i: v_{p_i}(disc K)} by the conductor-discriminant formula (exact integers)."""
    N = degree(F); out = {}
    for p, n in F:
        s = 0
        for a in range(n + 1):
            nch = 1 if a == 0 else p ** a - p ** (a - 1)            # characters of B_{p,n} of order p^a
            ce = 0 if a == 0 else (a + 2 if p == 2 else a + 1)      # their conductor exponent
            s += nch * ce
        out[p] = s * (N // p ** n)
    return out

def group_H(F):
    m = conductor(F)
    return [a for a in range(1, m) if math.gcd(a, m) == 1 and in_H(a, F)]

def coset_reps(F):
    m = conductor(F); H = group_H(F); seen = set(); reps = []
    for a in range(1, m):
        if math.gcd(a, m) != 1 or a in seen: continue
        reps.append(a)
        for h in H: seen.add(a * h % m)
    return reps

class Cyclo:
    """Z[zeta_m] = Z[x]/Phi_m(x), exact."""
    def __init__(self, m):
        self.m = m; self.Phi = fmpz_poly.cyclotomic(m)
    def red(self, a): return a % self.Phi
    def mul(self, a, b): return (a * b) % self.Phi
    def monomial_sum(self, exps):
        c = [0] * self.m
        for e in exps: c[e % self.m] += 1
        return fmpz_poly(c) % self.Phi

def gaussian_period(F):
    """eta = sum_{h in H} zeta_m^h as an element of Z[zeta_m] (exact), plus the cyclotomic ring."""
    m = conductor(F); Z = Cyclo(m)
    return Z, Z.monomial_sum(group_H(F))

def period_minpoly(F, prec=None):
    """Minimal polynomial T0 of eta: numerically from the conjugates eta_a = sum_h cos(2 pi a h/m)*... (arb),
       rounded, then T0(eta)=0 checked EXACTLY in Z[zeta_m]; the conjugates are checked (rigorously) to be pairwise
       distinct, so deg minpoly(eta) = N and T0 (monic, degree N) is the minimal polynomial."""
    m = conductor(F); N = degree(F); H = group_H(F); reps = coset_reps(F)
    assert len(reps) == N and len(H) * N == sum(1 for a in range(1, m) if math.gcd(a, m) == 1)
    assert (m - 1) in H                                             # -1 in H: K is totally real
    old = ctx.prec; ctx.prec = prec or (64 + 8 * N * max(1, int(math.log2(len(H) + 1)) + 2))
    tw = 2 * arb.pi() / m
    roots = [sum((tw * (a * h % m)).cos() for h in H) for a in reps]   # real since -1 in H
    for i in range(N):
        for j in range(i):
            assert not roots[i].overlaps(roots[j]), "conjugates not separated"
    P = arb_poly.from_roots(roots)
    co = []
    for c in P.coeffs():
        u = c.unique_fmpz(); assert u is not None, "coefficient not determined; raise precision"
        co.append(int(u))
    ctx.prec = old
    T0 = fmpz_poly(co)
    Z, eta = gaussian_period(F)
    assert horner_exact(Z, T0, eta) == 0
    return T0, Z, eta, roots

def horner_exact(Z, P, y, D=1):
    """D^deg P * P(y/D) in Z[zeta_m], for an integer polynomial P and y in Z[zeta_m] (exact)."""
    c = [int(x) for x in P.coeffs()]; d = len(c) - 1
    r = fmpz_poly([0]); Dp = 1
    for k in range(d, -1, -1):
        r = Z.mul(r, y) + fmpz_poly([c[k] * Dp]); Dp *= D
    # r = sum c_k y^k D^{d-k}
    return Z.red(r)

def eval_rational_at(Z, num, D, y):
    """D * (num(y)/D) = num(y) in Z[zeta_m] for an integer-coefficient poly num (list low->high)."""
    r = fmpz_poly([0])
    for a in reversed(num): r = Z.mul(r, y) + fmpz_poly([a])
    return r

def small_primes(n):
    if n < 2: return []
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def mult_order(l, N):
    k, x = 1, l % N
    while x != 1: x = x * l % N; k += 1
    return k

def hmin(F):
    """min over primes l not dividing N of l^{ord_N(l)} (Masley's rank corollary)."""
    N = degree(F); best = None
    for l in small_primes(10 ** 5):
        if best is not None and l >= best: break
        if N % l == 0: continue
        v = l ** mult_order(l, N)
        if best is None or v < best: best = v
    return best

FIELDS = [[(2,1),(3,1)],[(2,2),(3,1)],[(2,1),(5,1)],[(3,1),(5,1)],[(2,1),(7,1)],[(2,1),(3,2)],[(2,2),(5,1)],
          [(3,1),(7,1)],[(2,3),(3,1)],[(2,1),(11,1)],[(2,2),(7,1)],[(2,2),(3,2)],[(2,1),(13,1)]]
