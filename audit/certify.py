import sys,time
from fields import *
m=int(sys.argv[1])
P=Qm(m)
t=time.time()
bnf=pari.bnfinit(pari.nfinit([P,[p for p in factorint(m)]]),1)
print(m,'h(GRH)=',bnf.bnf_get_no(),flush=True)
r=pari.bnfcertify(bnf)
print(m,'bnfcertify ->',r,'time',time.time()-t,flush=True)
