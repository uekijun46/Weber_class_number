#!/usr/bin/env python3
# Audit k_41, k_43: tampered certificates for relations/span.py (written to tamper/).  The base is the header and the first 30 rows.
# After audit/k37/mk_tamper.py, with the rows of the present certificates.     python3 mk_tamper.py p
import sys
from fieldp import Field, read_cert, CERT, need
p=int(sys.argv[1]); F=Field(p); src=CERT(p)
lines=[l.rstrip("\n") for l in open(src)][:32]; hdr=lines[:2]; body=lines[2:]; cert=read_cert(src,p)
def w(name,ls): open("tamper/k%d_%s.txt"%(p,name),"w").write("\n".join(ls)+"\n")
r=body[2]; c,coef,fac=cert[2]; cs,co,fa=[t.strip() for t in r.split("|")]
big=[(q,e) for q,e in fac if e==1 and q!=p]; q1,q2=big[0][0],big[1][0]; i0,a0=coef[0]
other=next(q for q in sorted(set(q for _,_,fc in cert[:30] for q,e in fc)) if q not in dict(fac) and q!=p)
w("t00_untouched_30rows",hdr+body)
def rep(new): return hdr+body[:2]+[new]+body[3:]
w("t01_exponent",rep("%s | %s | %s"%(cs,co,fa.replace("%d^1"%q1,"%d^2"%q1,1))))
w("t02_coef_sign",rep("%s | %s | %s"%(cs,co.replace("%d:%d"%(i0,a0),"%d:%d"%(i0,-a0),1),fa)))
w("t03_c_plus1",rep("%d | %s | %s"%(c+1,co,fa)))
w("t04_composite",rep("%s | %s | %s"%(cs,co," ".join("%d^%d"%(q,e) for q,e in fac if q not in (q1,q2))+" %d^1"%(q1*q2))))
w("t05_header_other_p",["# field %d 1"%(84-p)]+hdr[1:]+body)
w("t06_header_layer2",["# field %d 2"%p]+hdr[1:]+body)
w("t07_nonsplit_2",rep("2 |  | 2^%d"%p))
w("t08_prime_swapped",rep("%s | %s | %s"%(cs,co,fa.replace("%d^1"%q1,"%d^1"%other,1))))
w("t09_extra_factor_p",rep(r+" %d^1"%p))
w("t10_duplicate_prime",rep("%s | %s | %s"%(cs,co,fa.replace("%d^1"%q1,"%d^1 %d^1"%(q1,q1),1)) if dict(fac)[q1]==1 else r))
w("t11_index_p",rep("%s | %s | %s"%(cs,co.replace("%d:%d"%(i0,a0),"%d:%d"%(i0+p,a0),1),fa)))
w("t12_shifted_indices",rep("%s | %s | %s"%(cs," ".join("%d:%d"%((i+1)%p,a) for i,a in coef),fa)))    # tau(alpha): a true relation (harmless, must pass)
# alpha = beta^2 for beta = the third row: a true element with a true factorization, but every prime ideal has exponent 2 or more
v=F.elem(c,coef); A=F.mul_matrix(c,coef); sq=[sum(v[k]*A[k][t] for k in range(p)) for t in range(p)]
w("t13_square_of_row3",rep("%d | %s | %s"%(sq[0]," ".join("%d:%d"%(i,sq[1+i]) for i in range(p-1) if sq[1+i])," ".join("%d^%d"%(q,2*e) for q,e in fac))))
# the same element with eta_i replaced by eta_(2i) (another numbering of the periods): the norm changes, must fail
w("t14_other_numbering",rep("%s | %s | %s"%(cs," ".join("%d:%d"%((2*i)%p,a) for i,a in coef),fa)))
print("k_%d: written; row 3 = %s"%(p,r))
