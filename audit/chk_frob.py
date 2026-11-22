from fields import *
from miller import *
from sympy import primerange
for m in [6,10,12,14,15,18,20,21]:
    P=Qm(m); nf=pari.nfinit([P,[p for p in factorint(m)]+[2,3,5,7]])
    bad=0
    S=dict()
    for q in primerange(2,200):
        dec=pari.idealprimedec(nf,q)
        fs=[(int(pr[2]),int(pr[3])) for pr in dec]  # e,f
        if q in factorint(m):
            v=factorint(m)[q]; mp_=m//q**v
            f=frob_order_lcm(q,mp_) if mp_>1 else 1
            exp=(q**v,f)
        else:
            exp=(1,frob_order_lcm(q,m))
        if any((e,f)!=exp for e,f in fs) or len(fs)*exp[0]*exp[1]!=m: bad+=1; print('BAD',m,q,fs,exp)
    print(m,'ok' if not bad else bad)
