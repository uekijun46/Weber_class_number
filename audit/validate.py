from a5 import *
from sympy import factorint
B3=Bp(3,1); B5=Bp(5,1)
res=[]; agree=0; dis=0
for p,Bpol in [(3,B3),(5,B5)]:
  for d in range(2,400):
    if any(e>1 for e in factorint(d).values()) or d%p==0: continue
    k=pari.bnfinit(pari(f'x^2-{d}'),1)
    if int(k.bnf_get_no())%p==0: continue
    r=fm(k,p); t=int(r[0]); rk=int(r[1])
    if t<2: continue
    pred = rk < t-1
    K=pari.polredbest(pari.polcompositum(pari(f'x^2-{d}'),Bpol)[0])
    bk=pari.bnfinit(K,1)
    act = int(bk.bnf_get_no())%p==0
    if pred==act: agree+=1
    else: dis+=1; print('DISAGREE',p,d,t,rk,bk.bnf_get_cyc())
    if pred: res.append((p,d,str(bk.bnf_get_cyc())))
print('agree',agree,'disagree',dis); print('predicted divisible cases:',res)
