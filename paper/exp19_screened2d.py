"""Screening on the 5-knob 2D study: does it move the 440x calibration cost?"""
import sys, os, itertools, pickle, numpy as np
from knobs import run_protocol
from knobs.problems import Burgers2D
b=Burgers2D(nu=0.05,T=1.0); EXACT=b.exact_observables()
tensor=pickle.load(open(os.path.join(os.path.dirname(__file__),"data","burgers2d_tensor.pkl"),"rb"))
LAD=[[16,32,64,128],[10,20,40,80],[1e-2,1e-5,1e-8],[1e-2,1e-5,1e-8],[1e-1,1e-2,1e-3]]
KINDS=["discretisation","discretisation","nested-solver","nested-solver","cost-only"]
HS=[[2*np.pi/n for n in LAD[0]],[1.0/s for s in LAD[1]]]+[list(l) for l in LAD[2:]]
ALL=list(itertools.product(*[range(len(l)) for l in LAD]))
print("obs        target    mode       safe    penalty(med/p90)  calib/user-run  runs  frozen")
for ob in ["energy","gradient"]:
    QS=EXACT[ob]
    look=lambda i,ob=ob:(tensor[tuple(i)][0][ob],tensor[tuple(i)][1])
    tr=lambda i:abs(tensor[i][0][ob]-QS)
    for er in [1e-2,3e-3]:
        eps=er*abs(QS)
        feas=[i for i in ALL if tr(i)<=eps]; c_opt=min(tensor[i][1] for i in feas)
        for use in [False,True]:
            safe=0;pens=[];cal=[];rn=[];fz=[];rep=0
            for a in ALL:
                r=run_protocol(look,LAD,a,eps,hs=HS,kinds=KINDS,use_screen=use)
                p=r["proposal"]; safe+=tr(p)<=eps
                pens.append(tensor[p][1]/c_opt); cal.append(r["protocol_cost"]/tensor[a][1])
                rn.append(r["runs"]); fz.append(len(r["frozen"])); rep+=r["repairs"]
            pens=np.array(pens)
            print(f"{ob:<10} {er:.0e}   {'screened' if use else 'baseline':<9} "
                  f"{100*safe/len(ALL):5.1f}%  {np.median(pens):5.2f}/{np.percentile(pens,90):5.2f}"
                  f"        {np.median(cal):8.1f}x  {np.median(rn):5.0f}   {np.median(fz):.1f}"
                  + (f"  repairs={rep}" if rep else ""))
