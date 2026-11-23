import glob, collections, sys
from sympy import primerange
files=sorted(glob.glob('/home/claude/weber_class_number/B72/cert*.txt'))
cnt=collections.Counter(); headers=collections.Counter(); bad=0; perfile={}
for fn in files:
    n=0
    for line in open(fn):
        if not line.strip(): continue
        if line.startswith('#'):
            headers[line.split()[0]]+=1; continue
        p=int(line.split(';',1)[0]); cnt[p]+=1; n+=1
    perfile[fn]=n
print(len(files),'files; rows',sum(cnt.values()),'distinct',len(cnt),'headers',headers)
dups={p:c for p,c in cnt.items() if c>1}; print('dup primes',len(dups), sorted(dups)[:10])
X=22*10**6
mine=[p for p in primerange(2,X+1) if pow(p,6,343)==1]
print('my count',len(mine),'max',mine[-1])
S=set(cnt); M=set(mine)
print('in certs not mine',sorted(S-M)[:20],'missing',len(M-S),sorted(M-S)[:20])
open('cert_primes.txt','w').write('\n'.join(map(str,sorted(S))))
open('my_primes.txt','w').write('\n'.join(map(str,mine)))
