"""Can a knob's structural class be inferred from its own cheap sweep?

If yes, the 'generic' tool needs only a parameter list and ladders.  If no, the
user must declare the class and the guarantee becomes conditional on that
declaration.  Ground-truth classes are known here, so this is checkable.

Signature features, all from one OAT sweep at one anchor:
  err_range     max |dQ| over the ladder, relative to the anchor value
  plateau       fraction of the FINE half where refining changes nothing
  monotone      fraction of adjacent steps moving the observable one way
  cost_slope    d log(cost) / d (ladder index)
  cost_interior 1 if cost has an interior minimum (neither end is cheapest)
"""
import pickle, os, time, numpy as np
from errbudget.problems import Burgers2D

NAMES = ["N","steps","newton_tol","gmres_tol","ilu_drop","gmres_restart",
         "precond_freq","newton_max"]
LAD  = [[16,32,64,128],[10,20,40,80],[1e-2,1e-5,1e-8],[1e-2,1e-5,1e-8],
        [1e-1,1e-2,1e-3],[10,30,100],[20,5,1],[2,5,20]]
TRUE = ["discretisation","discretisation","nested-solver","nested-solver",
        "cost-only","nested-solver","cost-only","nested-solver"]

b = Burgers2D(nu=0.05, T=1.0)
CACHE=os.path.join(os.path.dirname(__file__),"data","classify.pkl")
store = pickle.load(open(CACHE,"rb")) if os.path.exists(CACHE) else {}
def get(idx):
    idx=tuple(idx)
    if idx not in store:
        th={NAMES[i]:LAD[i][j] for i,j in enumerate(idx)}
        obs,C,w,_=b.run(th); store[idx]=(obs,C)
    return store[idx]

def signature(anchor, ob):
    base, _ = get(anchor)
    rows=[]
    for k in range(8):
        Q=[];Cc=[]
        for j in range(len(LAD[k])):
            cc=tuple(j if i==k else anchor[i] for i in range(8))
            o,C=get(cc); Q.append(o[ob]); Cc.append(C)
        Q=np.array(Q); Cc=np.array(Cc)
        d=np.abs(Q-Q[-1])
        rng=d.max()/ (abs(Q[-1])+1e-300)
        half=len(Q)//2
        fine=d[half:]
        plateau=float(np.mean(fine <= 1e-3*max(d.max(),1e-300)))
        diffs=np.diff(Q)
        mono=float(np.mean(np.sign(diffs)==np.sign(diffs[np.argmax(np.abs(diffs))]))) \
             if np.any(diffs) else 1.0
        lc=np.log(Cc); slope=float(np.polyfit(np.arange(len(lc)),lc,1)[0])
        interior=int(0 < int(np.argmin(Cc)) < len(Cc)-1)
        rows.append(dict(knob=NAMES[k], rng=rng, plateau=plateau, mono=mono,
                         slope=slope, interior=interior))
    return rows

def classify(r, tau_err=1e-6):
    """Rules stated in advance, not fitted to the answer."""
    if r["rng"] < tau_err:
        return "cost-only"
    if r["plateau"] >= 0.5:
        return "nested-solver"      # error vanishes before the ladder ends
    return "discretisation"         # keeps paying off as you refine

t0=time.time()
for ob in ["energy","gradient"]:
    for anchor,lab in [((1,1,1,1,1,1,1,1),"mid"),((2,1,0,1,1,0,1,1),"off-centre")]:
        rows=signature(anchor,ob)
        print(f"\n=== {ob}, anchor {lab} ===")
        print(f"  {'knob':<15}{'err_range':>11}{'plateau':>9}{'mono':>7}"
              f"{'cost_slope':>12}{'interior':>9}   {'inferred':<15}{'true':<15}")
        ok=0
        for r,t in zip(rows,TRUE):
            c=classify(r); ok+= (c==t)
            flag="" if c==t else "   <-- MISMATCH"
            print(f"  {r['knob']:<15}{r['rng']:>11.2e}{r['plateau']:>9.2f}{r['mono']:>7.2f}"
                  f"{r['slope']:>12.3f}{r['interior']:>9d}   {c:<15}{t:<15}{flag}")
        print(f"  correct: {ok}/8")
print(f"\n{len(store)} runs, {time.time()-t0:.0f}s")
pickle.dump(store, open(CACHE,"wb"))
