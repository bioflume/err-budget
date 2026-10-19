"""Does the screening exemption rule hold up with eight knobs?

Prediction from the non-transferability law: the knobs whose influence fails to
transfer from a cheap screening configuration to a production one are exactly the
NESTED ITERATIVE-SOLVER knobs -- newton_tol, newton_max, gmres_tol, gmres_restart
-- because changing the discretisation changes the problem they are solving.
Discretisation knobs (N, steps) and pure-cost knobs (ilu_drop, precond_freq)
should transfer.  Falsifiable either way.
"""
import sys, os, time, pickle, os, numpy as np
from knobs.problems import Burgers2D

NAMES = ["N","steps","newton_tol","gmres_tol","ilu_drop","gmres_restart",
         "precond_freq","newton_max"]
LAD  = [[16,32,64,128],[10,20,40,80],[1e-2,1e-5,1e-8],[1e-2,1e-5,1e-8],
        [1e-1,1e-2,1e-3],[10,30,100],[20,5,1],[2,5,20]]
KIND = ["discretisation","discretisation","nested solver","nested solver",
        "cost-only","nested solver","cost-only","nested solver"]

b = Burgers2D(nu=0.05, T=1.0)
EXACT = b.exact_observables()
CACHE = os.path.join(os.path.dirname(__file__),"data","screen8.pkl")
store = pickle.load(open(CACHE,"rb")) if os.path.exists(CACHE) else {}

def get(idx):
    idx = tuple(idx)
    if idx not in store:
        th = {NAMES[i]: LAD[i][j] for i,j in enumerate(idx)}
        obs, C, w, _ = b.run(th)
        store[idx] = (obs, C)
    return store[idx]

rng = np.random.default_rng(7)
COMBOS = [tuple(int(rng.integers(0,2)) for _ in range(6)) for _ in range(24)]
CHEAP, EXP = (0,0), (2,2)

def influence(anchor):
    """Normalised |change in observable| when each knob is refined one step."""
    base = get(anchor)[0]
    out = {}
    for ob in ["energy","gradient"]:
        s=[]
        for k in range(8):
            j = min(anchor[k]+1, len(LAD[k])-1)
            if j == anchor[k]: s.append(0.0); continue
            cc = tuple(j if i==k else anchor[i] for i in range(8))
            s.append(abs(get(cc)[0][ob]-base[ob]))
        s=np.array(s); t=s.sum()
        out[ob] = s/t if t>0 else s
    return out

t0=time.time()
res = {ob: {"fn":np.zeros(8), "fp":np.zeros(8)} for ob in ["energy","gradient"]}
n=0
for ci,c in enumerate(COMBOS):
    a_cheap = tuple(list(CHEAP)+list(c))
    a_exp   = tuple(list(EXP)+list(c))
    ic, ie = influence(a_cheap), influence(a_exp)
    n+=1
    for ob in res:
        for k in range(8):
            if ic[ob][k]<0.05 and ie[ob][k]>=0.05: res[ob]["fn"][k]+=1
            if ic[ob][k]>=0.05 and ie[ob][k]<0.05: res[ob]["fp"][k]+=1
    if ci%6==0:
        print(f"  {ci+1}/{len(COMBOS)} combos, {len(store)} runs, {time.time()-t0:.0f}s", flush=True)
pickle.dump(store, open(CACHE,"wb"))

print(f"\n{'='*80}\n8-knob screening: cheap (N=16,steps=10) -> production (N=64,steps=40)")
print(f"{n} parameter combinations, {len(store)} solver runs, {time.time()-t0:.0f}s\n")
for ob in ["energy","gradient"]:
    print(f"  {ob}:")
    print(f"    {'knob':<15} {'class':<16} {'missed by screen':>17} {'kept needlessly':>16}")
    for k,nm in enumerate(NAMES):
        print(f"    {nm:<15} {KIND[k]:<16} {int(res[ob]['fn'][k]):>10}/{n:<6} "
              f"{int(res[ob]['fp'][k]):>9}/{n}")
    fn=res[ob]["fn"]
    nested=[k for k in range(8) if KIND[k]=="nested solver"]
    other=[k for k in range(8) if KIND[k]!="nested solver"]
    print(f"    -> misses among nested-solver knobs: {int(fn[nested].sum())}"
          f"   among all others: {int(fn[other].sum())}")

# cost of screening vs production
cc = get(tuple(list(CHEAP)+list(COMBOS[0])))[1]
ce = get(tuple(list(EXP)+list(COMBOS[0])))[1]
print(f"\n  a screening run costs {cc:.3g}, a production run {ce:.3g} -> {ce/cc:.0f}x cheaper")
