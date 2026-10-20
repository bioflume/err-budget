"""Consolidated headline numbers under the final configuration."""
import sys, os, itertools, pickle, numpy as np
from errbudget import run_protocol
from errbudget.problems import Burgers1D
from errbudget.problems import Burgers2D

def study(name, exact, tensor, LAD, HS, obs_list, targets, rel=False):
    ALL=list(itertools.product(*[range(len(l)) for l in LAD]))
    print(f"\n{'='*92}\n{name}   ({len(ALL)} starting points)")
    print("obs          target      safe     penalty(med/p90/max)   already-ok: saved   "
          "missed: fixed   calib")
    for ob in obs_list:
        QS=exact[ob]
        look=lambda i,ob=ob:(tensor[tuple(i)][0][ob],tensor[tuple(i)][1])
        tr=lambda i:abs(tensor[i][0][ob]-QS)
        for t in targets:
            eps = t*abs(QS) if rel else t
            feas=[i for i in ALL if tr(i)<=eps]
            if not feas: continue
            c_opt=min(tensor[i][1] for i in feas)
            safe=0; pens=[]; ovh=[]; saved=[]; nbad=0; fixed=0
            for a in ALL:
                r=run_protocol(look,LAD,a,eps,hs=HS)
                p=r["proposal"]; ok=tr(p)<=eps; safe+=ok
                pens.append(tensor[p][1]/c_opt); ovh.append(r["protocol_cost"]/tensor[a][1])
                if tr(a)<=eps: saved.append(tensor[a][1]/tensor[p][1])
                else: nbad+=1; fixed+=ok
            pens=np.array(pens); sv=np.array(saved) if saved else np.array([np.nan])
            lab=f"{t:.0e}" + ("rel" if rel else "")
            print(f"{ob:<12} {lab:<10} {100*safe/len(ALL):5.1f}%  "
                  f"{np.median(pens):5.2f}/{np.percentile(pens,90):5.2f}/{pens.max():6.2f}   "
                  f"{np.nanmedian(sv):6.2f}x ({len(saved):3d})   "
                  f"{fixed:3d}/{nbad:<3d}       {np.median(ovh):6.1f}x")

b1=Burgers1D(nu=0.05,T=1.0)
t1=pickle.load(open(os.path.join(os.path.dirname(__file__),"data","burgers_tensor.pkl"),"rb"))
Ns,Ss=[32,64,128,256,512],[25,50,100,200,400]; NT=GT=[1e-2,1e-4,1e-6,1e-8,1e-10]
study("L2  1D viscous Burgers, Cole-Hopf truth, 4 knobs", b1.exact_observables(), t1,
      [Ns,Ss,NT,GT],
      [[2*np.pi/n for n in Ns],[1.0/s for s in Ss],list(NT),list(GT)],
      ["energy","dissipation"], [1e-2,3e-3,1e-3])

b2=Burgers2D(nu=0.05,T=1.0)
t2=pickle.load(open(os.path.join(os.path.dirname(__file__),"data","burgers2d_tensor.pkl"),"rb"))
Ns,Ss=[16,32,64,128],[10,20,40,80]; NT=GT=[1e-2,1e-5,1e-8]; DR=[1e-1,1e-2,1e-3]
study("L3  2D Burgers, manufactured solution, 5 knobs", b2.exact_observables(), t2,
      [Ns,Ss,NT,GT,DR],
      [[2*np.pi/n for n in Ns],[1.0/s for s in Ss],list(NT),list(GT),list(DR)],
      ["energy","gradient"], [1e-2,3e-3], rel=True)
