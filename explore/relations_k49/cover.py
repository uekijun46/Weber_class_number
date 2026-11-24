import pickle, sys, time
from proto import split, peel, work
from multiprocessing import Pool
rows=pickle.load(open("rows.pkl","rb"))
if len(sys.argv)>1:
    n=int(sys.argv[1]); s0=int(sys.argv[2]); t0=time.time()
    with Pool(2) as pool: res=pool.map(work,[(s0,n//2),(s0+1,n//2)])
    for out,tt,tf in res: rows+=out
    pickle.dump(rows,open("rows.pkl","wb")); print("added %d elements in %.0f s"%(n,time.time()-t0))
seed=(0,(0,),(1,),((19,1),))
t1=time.time(); cert,order=peel([seed]+rows)
X=22*10**6
sv=bytearray([1])*(X+1)
for i in range(2,int(X**.5)+1):
    if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
sp=[q for q in range(2,X+1) if sv[q] and split(q)]
print("rows: %d; certified primes in total: %d (peeling %.0f s)"%(len(rows),len(cert),time.time()-t1))
for b in (10**5,10**6,5*10**6,22*10**6):
    print("   split primes <= %d: %d of %d certified"%(b,sum(1 for q in sp if q<=b and q in cert),sum(1 for q in sp if q<=b)))
