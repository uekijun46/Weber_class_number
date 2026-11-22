import sys,time
from fields import *
m=int(sys.argv[1])
P=Qm(m)
t=time.time()
nf=pari.nfinit([P,[p for p in factorint(m)]])
bnf=pari.bnfinit(nf,1)
print(m,'deg',pari.poldegree(P),'clgp',bnf.bnf_get_cyc(),'time',time.time()-t,flush=True)
pari.writebin(f'bnf{m}.bin',bnf)
