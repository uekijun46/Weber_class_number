# Checks the results of the second check of the range files (workflow recheck.yml).
#   python3 check_recheck.py DIR          (DIR = relations/k3_4, with the subdirectory recheck/)
# For every committed checksum file rows_A_B.sha256 there must be a result recheck/recheck_A_B.json in which
#   - the independent check (audit/relations_recheck.py) and verify.py both passed for the field and the range [A,B),
#   - the checksum of the checked file equals the committed checksum,
#   - the number of split primes agrees between the independent check, verify.py and summary_A_B.json.
# It also reports which versions of the two programs (sha256) and which runs produced the results.
import sys, json, glob, os
d=sys.argv[1]; bad=[]; vers=set(); n=0; rows=0; primes=0; field=None
def need(cond,msg):
    if not cond: bad.append(msg)
for f in sorted(glob.glob(os.path.join(d,"rows_*.sha256"))):
    A,B=[int(x) for x in os.path.basename(f)[5:-7].split("_")]; n+=1; tag="[%d,%d)"%(A,B)
    sha=open(f).read().split()[0]
    try:
        r=json.load(open(os.path.join(d,"recheck","recheck_%d_%d.json"%(A,B)))); s=json.load(open(os.path.join(d,"summary_%d_%d.json"%(A,B))))
        a,b=r["recheck"],r["verify"]
        if field is None: field=s["field"]
        need(a["field"]==b["field"]==s["field"]==field,tag+": field")
        need(a["range"]==b["range"]==s["range"]==[A,B],tag+": range")
        need(a["rows_sha256"]==sha,tag+": checksum of the rows")
        need(type(s["split_primes"]) is int and a["split_primes"]==b["split_primes"]==s["split_primes"],tag+": number of split primes")
        need(a["rows"]==b["rows"]==s["rows"]>=s["split_primes"],tag+": number of rows")
        vers.add((r["run"],r["commit"],r["verify_sha256"],r["recheck_sha256"])); rows+=a["rows"]; primes+=a["split_primes"]
    except (OSError,ValueError,KeyError,TypeError) as e:
        bad.append(tag+": result missing or unreadable (%s)"%type(e).__name__)
extra=[f for f in glob.glob(os.path.join(d,"recheck","recheck_*.json")) if not os.path.exists(os.path.join(d,"rows_"+os.path.basename(f)[8:-5]+".sha256"))]
need(not extra,"results without a checksum file: %s"%extra[:3]); need(n>0,"no checksum files")
for v in sorted(vers): print("run %s, commit %s, verify.py sha256 %s, relations_recheck.py sha256 %s"%v)
for m in bad[:20]: print("BAD",m)
if bad: raise SystemExit("FAILED: %d problems"%len(bad))
print("RECHECK COMPLETE: field k_{%d^%d}, %d range files, %d rows, %d split primes; both checks passed for every file"%(field[0],field[1],n,rows,primes))
