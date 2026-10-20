"""Build the L2 tensor (data/burgers_tensor.pkl): 5^4 = 625 runs of the 1D Burgers solver."""
import itertools, pickle, time, os
from errbudget.problems import Burgers1D
b=Burgers1D(nu=0.05,T=1.0)
LAD=[[32,64,128,256,512],[25,50,100,200,400],[1e-2,1e-4,1e-6,1e-8,1e-10],[1e-2,1e-4,1e-6,1e-8,1e-10]]
NAMES=["N","steps","newton_tol","gmres_tol"]
t={}; t0=time.time()
for idx in itertools.product(*[range(len(l)) for l in LAD]):
    obs,C,w,_=b.run({NAMES[j]:LAD[j][v] for j,v in enumerate(idx)}); t[idx]=(obs,C,w)
pickle.dump(t,open(os.path.join(os.path.dirname(__file__),"data","burgers_tensor.pkl"),"wb")); print("L2 built",len(t),time.time()-t0)
