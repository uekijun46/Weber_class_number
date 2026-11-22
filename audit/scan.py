from a5 import *
from sympy import primerange
Ks={'Q(sqrt2)':2,'B22':4,'B31':3,'B32':9,'B51':5,'B71':7,'B11,1':11,'B13,1':13,'Q(6)':6}
for name,mm in Ks.items():
    b=bnfK(mm)
    out=[]; cnt=0
    for p in primerange(3,400):
        if mm%p==0: continue
        r=fm(b,p); t=int(r[0]); rk=int(r[1])
        if t>=2:
            cnt+=1
            if rk<t-1: out.append((p,t,rk))
    print(name,'h=',b.bnf_get_no(),'#p with t>=2:',cnt,'deficient:',out,flush=True)
