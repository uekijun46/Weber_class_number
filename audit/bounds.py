import json,sys
from miller import *
mp.mp.dps=18
claimed={6:(0,0,2.04),12:(0,0,2.75),10:(0,0,5.54),15:(5,449,5.88),14:(5,313,3.61),18:(0,0,4.68),20:(5,593,32.18),21:(5,863,29.83),24:(5,991,7.73),22:(141,29201,18.09),28:(60,13457,22.44),36:(80,27919,27.44),26:(1052,319729,21.67)}
ms=[int(a) for a in sys.argv[1:]]
for m in ms:
    k,qc,bc=claimed[m]
    try: split=json.load(open(f'split{m}.json'))
    except: split=[]
    use=[q for q in split if q<=qc]
    Smax = max(qc*2, 20000)
    base=principal_primes(m, max(Smax, 400000 if m==26 else Smax))
    cs=[1.5*i for i in range(1,30)] if m!=26 else [26,28,30,32,34,36,38,40]
    def bb(sp):
        best=(mp.inf,None)
        for c in cs:
            b,B=bound(m,c,base+[(q,m) for q in sp])
            if b<best[0]: best=(b,c)
        return best
    b,c=bb(use)
    print(m,'hmin',hmin(m),'claimed',(k,qc,bc),'mine with same split set: n=',len(use),'bound',float(b),'c',c,flush=True)
