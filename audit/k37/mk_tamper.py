#!/usr/bin/env python3
# Audit k_37: tampered certificates for relations/span.py (written to tamper/).  The base is the header and the first 40 rows.
import sys
from field37 import Field, read_cert
F=Field(37); src="../../relations/k37_1/span_cert.txt"
lines=[l.rstrip("\n") for l in open(src)][:42]; hdr=lines[:2]; body=lines[2:]
def w(name,ls): open("tamper/%s.txt"%name,"w").write("\n".join(ls)+"\n")
w("t00_untouched_40rows",hdr+body)
r=body[2]
if "1069^2" not in r: raise SystemExit("unexpected third row")
w("t01_exponent",hdr+body[:2]+[r.replace("1069^2","1069^1")]+body[3:])
w("t02_coef_sign",hdr+body[:2]+[r.replace("2:-1","2:1")]+body[3:])
w("t03_c_plus1",hdr+body[:2]+[r.replace("-1 |","0 |",1)]+body[3:])
w("t04_composite",hdr+body[:2]+[r.replace("1669^1 1723^1","%d^1"%(1669*1723))]+body[3:])
w("t05_header_41",["# field 41 1"]+hdr[1:]+body)
w("t06_header_37_2",["# field 37 2"]+hdr[1:]+body)
w("t07_nonsplit_2",hdr+body[:2]+["2 |  | 2^37"]+body[3:])
w("t08_prime_swapped",hdr+body[:2]+[r.replace("1723^1","1951^1")]+body[3:])
w("t09_extra_factor_37",hdr+body[:2]+[r+" 37^1"]+body[3:])
w("t10_duplicate_prime",hdr+body[:2]+[r.replace("1069^2","1069^1 1069^1")]+body[3:])
w("t11_index_37",hdr+body[:2]+[r.replace("22:1","37:1")]+body[3:])
w("t12_zero_coef",hdr+body[:2]+[r.replace("22:1","22:1 5:0")]+body[3:])
# alpha = beta^2 for beta = the third row: every prime ideal that contains alpha has exponent 2 or 4 (a true element with a true factorization)
c,coef,fac=read_cert(src)[2]; v=F.elem(c,coef); M=F.mul_matrix(v); sq=[sum(v[k]*M[k][t] for k in range(37)) for t in range(37)]
row="%d | %s | %s"%(sq[0]," ".join("%d:%d"%(i,sq[1+i]) for i in range(36) if sq[1+i])," ".join("%d^%d"%(r_,2*e) for r_,e in fac))
w("t13_square_of_row3",hdr+body[:2]+[row]+body[3:])
# alpha = 691 (norm 691^37, contained in all 37 prime ideals above 691: the entry is 0 in A; harmless)
w("t14_alpha_691_harmless",hdr+body+["691 |  | 691^37"])
print("written")
