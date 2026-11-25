import pickle, time, sys
import cypari2
pari=cypari2.Pari(); pari.allocatemem(10**9)
st=pickle.load(open("out/elim_state.pkl","rb"))
rows=st["rows"]; alive=st["alive"]; cols=st["cols"]
print(len(alive),len(cols),st["stuck"])
Phi=pari.polcyclo(37)
t0=time.time(); cnt=0; units=0; hist={}; seen={}
for i in alive[:400]:
    for q,v in rows[i].items():
        key=tuple(v)
        if key not in seen:
            seen[key]=abs(int(pari.polresultant(pari.Pol(list(reversed(v))),Phi)))
        nm=seen[key]; cnt+=1
        if nm==1: units+=1
        mx=max(abs(c) for c in v); hist[mx]=hist.get(mx,0)+1
print("entries",cnt,"units",units,"time",time.time()-t0, "distinct", len(seen))
print(sorted(hist.items())[:20])
import collections
print(collections.Counter(len(rows[i]) for i in alive).most_common(8))
