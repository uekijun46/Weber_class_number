#!/usr/bin/env python3
# Audit k_47: the least possible class number > 1.  If a prime l != 47 divides h(k_47), then l^f divides h (f = order of l modulo 47;
# rank theorem, the only proper subfield is Q); 47 does not divide h (Iwasawa).  So h > 1 implies h >= min over primes l != 47 of l^f.
def is_prime(x):
    if x<2: return False
    d=2
    while d*d<=x:
        if x%d==0: return False
        d+=1
    return True
p=47
def order(l):
    f=1
    while pow(l,f,p)!=1: f+=1
    return f
cand=[k*p+1 for k in range(1,8)]
print("numbers k*47+1, k=1..7:",[(c,"prime" if is_prime(c) else "composite") for c in cand])
first=next(c for c in (k*p+1 for k in range(1,100)) if is_prime(c)); print("least prime = 1 mod 47:",first)
vals=sorted((l**order(l),l,order(l)) for l in range(2,5000) if is_prime(l) and l!=p)
print("three smallest values of l^f for primes l < 5000:",vals[:3])
print("primes l < %d with their orders modulo 47:"%first," ".join("%d:%d"%(l,order(l)) for l in range(2,first) if is_prime(l) and l!=p))
bad=[(l,order(l)) for l in range(2,first) if is_prime(l) and l!=p and l**order(l)<first]
print("primes l < %d with l^f < %d: %s"%(first,first,bad))
print("l^2 >= %d for l >= %d, so for f >= 2 only l <= %d must be looked at; their orders: %s"%(first,next(l for l in range(2,100) if l*l>=first),next(l for l in range(2,100) if l*l>=first)-1,[(l,order(l)) for l in range(2,17) if is_prime(l)]))
print("multiplicative order must divide 46: possible f in {1,2,23,46}; primes of order 2 (l = -1 mod 47) below 5000:",[l for l in range(2,5000) if is_prime(l) and l%p==p-1][:5])
