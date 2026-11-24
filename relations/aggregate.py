# Combine the verified summaries of all ranges and compute the class number bound.
#   python3 aggregate.py DIR        (DIR contains summary_A_B.json files written from the output of verify.py)
# Each summary states: every split prime q with A <= q < B has a verified row, and gives rigorous lower bounds of the
# partial sums over these primes.  The ranges must be contiguous and start at 0.
import sys, json, glob, os
from flint import arb, ctx
from verify import K, bound, up, down, need
d=sys.argv[1]; S=[json.load(open(f)) for f in glob.glob(os.path.join(d,"summary_*.json"))]
need(S,"no summaries"); S.sort(key=lambda s:s["range"][0])
p,n=S[0]["field"]; F=K(p,n); end=0
for s in S:
    need(s["field"]==[p,n] and s["range"][0]==end,"ranges are not contiguous at %d"%end); end=s["range"][1]
tot=sum(s["split_primes"] for s in S); rows=sum(s["rows"] for s in S)
print("field k_{%d^%d}; %d ranges covering [0, %d); split primes with a verified row: %d; rows: %d"%(p,n,len(S),end,tot,rows))
logrd,ex=F.logrd(); hm=F.hmin(); ctx.prec=200; pi=arb.pi(); best=None
for c in sorted(int(c) for c in S[0]["sums"]):
    if not all(str(c) in s["sums"] for s in S): continue
    Ssum=arb(0)
    for s in S: Ssum+=arb(str(s["sums"][str(c)]))-arb("1e-12")           # decimal lower bounds, with a margin for the conversion
    I_up,_,_,_=bound(F,[],c)
    B=(pi/2+arb.const_euler()+(8*pi).log()-logrd-I_up+Ssum.lower()).lower()
    hb=(2*arb(c)*pi.sqrt()/(F.N*B)).upper() if B>0 else None
    print("c=%d: integral <= %.10f, prime sum >= %.10f, B >= %.7f"%(c,up(I_up,10),down(Ssum.lower(),10),down(B,7))+("" if hb is None else ", h < %.4f"%up(hb,4)))
    if hb is not None and (best is None or hb<best): best=hb
print("least possible class number > 1: %d"%hm)
print("VERDICT: PASS  (h < %.4f < %d)"%(up(best,4),hm) if best is not None and best<hm else "VERDICT: bound not below %d yet"%hm)
