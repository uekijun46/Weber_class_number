# Independent re-check of the B_{7,2} certificates with PARI (the main verifier B72/verify_cert.py uses FLINT and Z[zeta_343]).
# Each row  p ; D ; a_0 ... a_48  encodes  alpha = (1/D) sum a_i theta^i,  theta a root of T (header line "#T ...").
# Checks: T = polsubcyclo(343,49), monic, disc of the field = 7^138;  p prime and p^6 = 1 mod 343 (p splits completely);
#         alpha lies in the maximal order (nfinit computes it without GRH);  |Res(T, D*alpha)| = p * D^49, i.e. |N(alpha)| = p.
# Usage: python3 audit/b72_indep_check.py PART NPARTS      (run from the repository root)
import cypari2, sys, glob, time
part, nparts = int(sys.argv[1]), int(sys.argv[2])
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9, silent=True)
T0 = pari("polsubcyclo(343,49)"); assert pari.polcoef(T0, 49) == 1 and pari.nfdisc(T0) == pari(7)**138
nf = pari.nfinit(T0); t0 = time.time(); n = 0; k = 0; primes = set()
for fn in sorted(glob.glob("B72/cert*.txt")):
    T = None
    for line in open(fn):
        if line.startswith("#T"):
            T = pari.Polrev([int(x) for x in line.split()[1:]]); assert T == T0, fn
            continue
        if not line.strip() or line.startswith("#"): continue
        k += 1
        if k % nparts != part: continue
        assert T is not None
        ps, Ds, cs = line.split(";"); p = int(ps); D = int(Ds); A = pari.Polrev([int(x) for x in cs.split()])
        assert D >= 1 and pari.isprime(p) and pow(p, 6, 343) == 1, ("prime", p)
        assert pari.denominator(pari.nfalgtobasis(nf, A / D)) == 1, ("not integral", p)
        assert abs(int(pari.polresultant(T, A))) == p * D**49, ("norm", p)
        primes.add(p); n += 1
print("part %d/%d: %d rows checked, %d distinct primes, max %d, %.0fs: ALL OK" % (part, nparts, n, len(primes), max(primes), time.time() - t0))
