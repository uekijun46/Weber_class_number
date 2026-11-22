# Generate principality certificates for split primes of B_{7,2} (degree 49, conductor 343).
# Usage: python3 gen_cert.py LO HI OUTFILE [TIME_LIMIT_SEC] [WORKERS] [INDEX]
# With WORKERS and INDEX, only every WORKERS-th split prime is handled (for parallel runs).
# Primes already present in any cert*.txt of the current directory are skipped, so a run can be resumed.
# Each output row:  p ; D ; a_0 ... a_48   meaning alpha = (1/D) * sum a_i theta^i,  theta a root of T = polsubcyclo(343,49).
# bnfinit assumes GRH internally, but each row is checked afterwards by verify_cert.py without any hypothesis.
import cypari2, sys, time, os, glob
lo,hi,out=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
limit=float(sys.argv[4]) if len(sys.argv)>4 else 1e18
W=int(sys.argv[5]) if len(sys.argv)>6 else 1
R=int(sys.argv[6]) if len(sys.argv)>6 else 0
pari=cypari2.Pari(); pari.allocatemem(3*10**9, silent=True)
pari("default(realprecision,250)")   # fewer slow primes than at the default precision
t0=time.time()
if os.path.exists("bnf_B72.bin"):
    pari('b=read("bnf_B72.bin")')
else:
    pari("b=bnfinit(polsubcyclo(343,49),1)"); pari('writebin("bnf_B72.bin",b)')
print("bnf ready %.1fs, h(GRH)=%s"%(time.time()-t0,pari("b.no")),flush=True)
pari("""cert(p)={my(P=idealprimedec(b,p)[1],v=bnfisprincipal(b,P),g,d); if(#v[1]>0 && v[1]!=0, return(0)); g=lift(nfbasistoalg(b,v[2])); d=denominator(content(g)); [d, Vecrev(d*g,49)]}""")
done=set()
for fn in set(glob.glob("cert*.txt"))|({out} if os.path.exists(out) else set()):
    done|={int(l.split(";")[0]) for l in open(fn) if l.strip() and not l.startswith("#")}
new=not os.path.exists(out)
f=open(out,"a"); p=lo-1; n=0; tmax=0; k=0
if new: f.write("#T "+" ".join(str(x) for x in pari("Vecrev(b.pol)"))+"\n")
while time.time()-t0<limit:
    p=int(pari(f"nextprime({p+1})"))
    if p>hi: break
    if pow(p,6,343)!=1: continue
    k+=1
    if k%W!=R or p in done: continue
    t=time.time(); r=pari(f"cert({p})"); dt=time.time()-t; tmax=max(tmax,dt)
    if r==0: print("NON-PRINCIPAL?",p,flush=True); continue
    f.write(f"{p} ; {r[0]} ; "+" ".join(str(x) for x in r[1])+"\n"); f.flush(); n+=1
print("rows written:",n,"last prime tried:",p,"max time/prime %.1fs"%tmax,"total %.1fs"%(time.time()-t0),flush=True)
