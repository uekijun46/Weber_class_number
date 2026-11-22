import time
from fields import *
P=Qm(33); t=time.time()
bnf=pari.bnfinit(pari.nfinit([P,[3,11]]),1); print('Q(33) cyc',bnf.bnf_get_cyc(),time.time()-t,flush=True)
