import time
from fields import *
pari.allocatemem(4*10**8); pari.default("parisizemax",5*10**9)
P=Qm(62); print('deg',pari.poldegree(P),flush=True)
t=time.time()
nf=pari.nfinit([P,[2,31]]); print('nf',time.time()-t,flush=True)
bnf=pari.bnfinit(nf,1,[0.05,0.05]); print('cyc',bnf.bnf_get_cyc(),time.time()-t,flush=True)
