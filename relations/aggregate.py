# Final step for a field covered by ranges: combine the ranges and compute the bound of the class number.
#   python3 aggregate.py DIR [processes] [p n END]      (DIR contains the files summary_A_B.json written from the output of verify.py)
# With p n END, the program fails unless the field is k_{p^n} and the ranges cover exactly [0, END).
# The exit status is 0 only if the bound is below the least possible class number > 1.
#
# What is taken from a summary: only the statement "verify.py passed for the range [A,B)", that is, every split prime q with
# A <= q < B has a verified row.  The numbers in the summaries are not trusted.  This program computes again, from scratch,
#   - the list of the split primes in every range (with proven primality), and compares their number with the summary,
#   - rigorous lower bounds of the sums over these primes (ball arithmetic, results kept as integers),
# and it checks that the ranges are [A,B) with A < B, start at 0, and follow each other without a gap.
#
# The sum of Miller's lemma over a prime q is 2 * sum_m log(q) q^(-m/2) F(m log q) with F(x) = exp(-(x/c)^2)/cosh(x/2).
# Since cosh((m log q)/2) = (q^(m/2) + q^(-m/2))/2, the m-th term equals 4 log(q) exp(-(m log(q)/c)^2) / (q^m + 1).
# Every term is positive, so a sum over fewer values of m is a lower bound.
import sys, json, glob, os
from multiprocessing import Pool
from flint import arb, ctx
from verify import K, bound, up, down, need, split_primes, _frac
SCALE=10**30
def range_sums(job):
    p,n,A,B,cs=job; F=K(p,n); ctx.prec=128
    sp=split_primes(F,A,B); S={c:arb(0) for c in cs}; ic={c:1/arb(c)**2 for c in cs}
    for q in sp:
        L=arb(q).log(); L2=L*L; m=1; qm=q
        while True:
            w=4*L/(qm+1)
            for c in cs: S[c]+=w*(-(m*m)*L2*ic[c]).exp()
            if m==8 or qm>10**16: break
            m+=1; qm*=q
    low={}
    for c in cs:
        a,b=_frac(S[c].lower()); low[c]=a*SCALE//b                # floor(lower end * SCALE): an integer
    return A,B,len(sp),low
def main():
    d=sys.argv[1]; procs=int(sys.argv[2]) if len(sys.argv)>2 else (os.cpu_count() or 1)
    S=[json.load(open(f)) for f in glob.glob(os.path.join(d,"summary_*.json"))]
    need(S,"no summaries")
    for s in S:
        need(isinstance(s,dict) and all(k in s for k in ("field","range","rows","split_primes")),"summary: missing entries")
        need(len(s["range"])==2 and all(type(x) is int for x in s["range"]) and 0<=s["range"][0]<s["range"][1],"summary: the range must be [A,B) with 0 <= A < B")
        need(type(s["rows"]) is int and type(s["split_primes"]) is int and s["rows"]>=s["split_primes"]>=0,"summary: number of rows")
    S.sort(key=lambda s:s["range"][0])
    p,n=S[0]["field"]; F=K(p,n); end=0
    for s in S:
        need(s["field"]==[p,n],"summaries of different fields")
        need(s["range"][0]==end,"ranges are not contiguous at %d"%end); end=s["range"][1]
    if len(sys.argv)>3: need([p,n,end]==[int(x) for x in sys.argv[3:6]],"expected field and end: %s; found: k_{%d^%d}, [0, %d)"%(sys.argv[3:6],p,n,end))
    cs=[60,80,100,150,200]
    print("field k_{%d^%d}; %d ranges covering [0, %d)"%(p,n,len(S),end),flush=True)
    with Pool(procs) as pool:
        res={(A,B):(cnt,low) for A,B,cnt,low in pool.imap_unordered(range_sums,[(p,n,s["range"][0],s["range"][1],cs) for s in S])}
    tot=0; total={c:0 for c in cs}
    for s in S:
        A,B=s["range"]; cnt,low=res[(A,B)]
        need(cnt==s["split_primes"],"range [%d,%d): %d split primes found, the summary states %d"%(A,B,cnt,s["split_primes"]))
        for c in cs:
            total[c]+=low[c]
            if "sums" in s and str(c) in s["sums"]:                # consistency only; the recomputed value is the one that is used
                need(abs(arb(low[c])/SCALE-arb(str(s["sums"][str(c)])))<arb("1e-9"),"range [%d,%d): the sum in the summary disagrees with the recomputed sum"%(A,B))
        tot+=cnt
    print("split primes below %d, each with a verified row according to the summaries: %d (number recomputed here); rows: %d"%(end,tot,sum(s["rows"] for s in S)))
    print("prime sums: recomputed here from the list of split primes")
    logrd,ex=F.logrd(); hm=F.hmin(); ctx.prec=200; pi=arb.pi(); best=None
    for c in cs:
        S_low=(arb(total[c])/SCALE).lower()
        I_up,_,_,_=bound(F,[],c)
        B=(pi/2+arb.const_euler()+(8*pi).log()-logrd-I_up+S_low).lower()
        hb=(2*arb(c)*pi.sqrt()/(F.N*B)).upper() if B>0 else None
        print("c=%d: integral <= %s, prime sum >= %s, B >= %s"%(c,up(I_up,10),down(S_low,10),down(B,7))+("" if hb is None else ", h < %s"%up(hb,4)))
        if hb is not None and (best is None or hb<best): best=hb
    print("degree %d, discriminant %d^%d, least possible class number > 1: %d"%(F.N,p,ex,hm))
    if best is not None and best<hm: print("VERDICT: PASS  (h < %s < %d)"%(up(best,4),hm))
    else: raise SystemExit("VERDICT: bound not below %d yet"%hm)
if __name__=="__main__": main()
