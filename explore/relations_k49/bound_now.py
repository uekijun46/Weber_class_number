import pickle, sys, math
sys.path.insert(0,'/home/claude/weber_class_number/B72')
from proto import split, peel
import verify_cert as v
rows=pickle.load(open("rows.pkl","rb")); seed=(0,(0,),(1,),((19,1),))
cert,order=peel([seed]+rows)
ps=sorted(cert)
print("certified split primes:",len(ps),"; below 2.2e7:",sum(1 for q in ps if q<=22*10**6),"; below 1e9:",sum(1 for q in ps if q<=10**9),"; largest: %.1e"%ps[-1])
for cut in (22*10**6,10**8,10**9,10**10,None):
    S=[q for q in ps if cut is None or q<=cut]
    best=None
    for c in (60,80,100,120,150,200):
        I_up,S_low,B,hb=v.bound(S,c)
        if hb is not None and (best is None or float(hb)<best[0]): best=(float(hb),c,float(B))
    print("  S = certified primes <= %s (%d primes): best bound"%(cut,len(S)), ("h < %.2f at c=%d (B >= %.5f)"%best if best else "B <= 0"))
