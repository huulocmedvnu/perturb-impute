import numpy as np, pandas as pd, time
from scipy.stats import spearmanr
m=pd.read_csv("gwps.meta.tsv",sep="\t"); gn=pd.read_csv("gwps.genes.tsv",sep="\t")
Z=np.fromfile("gwps.z.f32",np.float32).reshape(3,len(m),len(gn)).astype(np.float64)
fo=pd.read_csv("folds.tsv",sep="\t").set_index("pert").fold.reindex(m.pert).values
te=np.where(fo==0)[0]; rf=np.where((fo!=0)&(fo>=0))[0]
t=time.time(); Y=Z[0][rf]; w,V=np.linalg.eigh(Y.T@Y); V=V[:,::-1]; print("eigh",time.time()-t)
def spa(v,tt):
    d=np.abs(tt)>4
    return spearmanr(np.abs(v[d]),np.abs(tt[d]))[0] if d.sum()>=5 else np.nan
H0,H1=Z[1][te],Z[2][te]
print("direct", np.nanmedian([spa(H0[i],H1[i]) for i in range(len(te))]))
for K in [20,160,500,1000,2000,4000]:
    P=H0@V[:,:K]@V[:,:K].T
    print("program K",K, round(np.nanmedian([spa(P[i],H1[i]) for i in range(len(te))]),3))
G=pd.read_csv("cvk.genes.tsv",sep="\t"); uset=set(G[(G.fold==0)&(G.panel=="designed_500")].gene)
u=np.array([j for j,g in enumerate(gn.gene) if g in uset]); print("common genes",len(u))
print("COMMON direct", round(np.nanmedian([spa(H0[i,u],H1[i,u]) for i in range(len(te))]),3))
for K in [20,160,500,2000]:
    P=(H0@V[:,:K])@V[u,:K].T
    print("COMMON program K",K, round(np.nanmedian([spa(P[i],H1[i,u]) for i in range(len(te))]),3))
k=pd.read_csv("cvk.knn.tsv",sep="\t"); k=k[(k.fold==0)&(k.panel=="designed_500")]; print("ridge fold0", k.ridge_spa.median())
