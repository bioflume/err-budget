"""Q1: can we screen out irrelevant knobs CHEAPLY and trust that at production cost?

If the set of knobs that matter can be identified at a coarse, cheap configuration
and still holds at an expensive one, then the number of knobs stops being the
obstacle: screen at low cost, then treat only the survivors properly.  This is a
transfer claim of exactly the kind that already failed once here, so it is tested,
not assumed.
"""
import sys, os, itertools, pickle, numpy as np

def shares(tensor, LAD, c, ob):
    """Normalised influence of each knob at config c: one step finer, each knob."""
    base = tensor[c][0][ob]
    s = []
    for k in range(len(LAD)):
        j = min(c[k]+1, len(LAD[k])-1)
        if j == c[k]:
            s.append(0.0); continue
        cc = tuple(j if i==k else c[i] for i in range(len(LAD)))
        s.append(abs(tensor[cc][0][ob]-base))
    s = np.array(s); t = s.sum()
    return s/t if t>0 else s

def screen_test(name, tensor, LAD, NAMES, obs_list, cheap_ix, exp_ix, tau=0.05):
    print(f"\n{'='*84}\n{name}   (threshold: a knob 'matters' if it holds >{tau:.0%} of influence)")
    nk=len(LAD)
    others=[range(len(LAD[i])) for i in range(2,nk)]
    for ob in obs_list:
        fn=np.zeros(nk); fp=np.zeros(nk); n=0
        for rest in itertools.product(*others):
            c=tuple(list(cheap_ix)+list(rest))
            C=tuple(list(exp_ix)+list(rest))
            sc, sC = shares(tensor,LAD,c,ob), shares(tensor,LAD,C,ob)
            n+=1
            for k in range(nk):
                if sc[k]<tau and sC[k]>=tau: fn[k]+=1     # missed at low cost: DANGEROUS
                if sc[k]>=tau and sC[k]<tau: fp[k]+=1     # kept needlessly: merely wasteful
        print(f"  {ob}:  {n} cheap/expensive pairs")
        print(f"    {'knob':<12} {'missed by cheap screen':>24} {'kept unnecessarily':>20}")
        for k,nm in enumerate(NAMES):
            print(f"    {nm:<12} {int(fn[k]):>14}/{n:<8} {int(fp[k]):>14}/{n}")

# ---- 2D
t2=pickle.load(open(os.path.join(os.path.dirname(__file__),"data","burgers2d_tensor.pkl"),"rb"))
LAD2=[[16,32,64,128],[10,20,40,80],[1e-2,1e-5,1e-8],[1e-2,1e-5,1e-8],[1e-1,1e-2,1e-3]]
screen_test("2D Burgers: screen at N=16,steps=10  ->  check at N=128,steps=80",
            t2, LAD2, ["N","steps","newton_tol","gmres_tol","ilu_drop"],
            ["energy","gradient"], (0,0), (2,2))

# ---- 1D
t1=pickle.load(open(os.path.join(os.path.dirname(__file__),"data","burgers_tensor.pkl"),"rb"))
LAD1=[[32,64,128,256,512],[25,50,100,200,400],
      [1e-2,1e-4,1e-6,1e-8,1e-10],[1e-2,1e-4,1e-6,1e-8,1e-10]]
screen_test("1D Burgers: screen at N=32,steps=25  ->  check at N=256,steps=200",
            t1, LAD1, ["N","steps","newton_tol","gmres_tol"],
            ["energy","dissipation"], (0,0), (3,3))

# ---- what does screening actually save?
print(f"\n{'='*84}\ncost of screening vs treating every knob at production resolution")
for nm,t,LAD,cheap,exp in [("2D",t2,LAD2,(0,0),(2,2)),("1D",t1,LAD1,(0,0),(3,3))]:
    nk=len(LAD)
    c=tuple(list(cheap)+[len(l)//2 for l in LAD[2:]])
    C=tuple(list(exp)+[len(l)//2 for l in LAD[2:]])
    print(f"  {nm}: one screening run costs {t[c][1]:.3g}, one production run {t[C][1]:.3g} "
          f"-> screening is {t[C][1]/t[c][1]:.0f}x cheaper per run")
