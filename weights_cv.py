import numpy as np, itertools
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
X=np.load('eps_inputs.npy'); y=np.load('eps_labels.npy')
vina=X[:,0]; dE=X[:,1]; f=X[:,2]; b=X[:,3]; sl=X[:,4]
dEo,fo,bo=-9.386,0.1281,484.7
RT=0.616033; EPS=0.02
def S(w1,w2,w3):
    s=w1*(dE-dEo) - w2*RT*np.log((f+EPS)/(fo+EPS)) - w3*RT*(b/sl - bo/sl)
    return np.where(dE>0.0, 1e6, s)
grid=[0.0,0.5,1.0,1.5,2.0,3.0]
combos=[c for c in itertools.product(grid,grid,grid) if c!=(0,0,0)]
sss=StratifiedShuffleSplit(n_splits=200, test_size=0.3, random_state=0)
test_aucs=[]; fixed_aucs=[]; chosen=[]
for tr,te in sss.split(X,y):
    best=None
    for c in combos:
        a=roc_auc_score(y[tr], -S(*c)[tr])
        if best is None or a>best[0]: best=(a,c)
    chosen.append(best[1])
    test_aucs.append(roc_auc_score(y[te], -S(*best[1])[te]))
    fixed_aucs.append(roc_auc_score(y[te], -S(1.0,1.0,1.5)[te]))
test_aucs=np.array(test_aucs); fixed_aucs=np.array(fixed_aucs)
print('splits=200 | test 30%% | grid %d combos' % len(combos))
print('train-chosen weights: test AUC mean %.3f [%.3f, %.3f]' % (test_aucs.mean(), np.percentile(test_aucs,2.5), np.percentile(test_aucs,97.5)))
print('fixed (1,1,1.5):     test AUC mean %.3f [%.3f, %.3f]' % (fixed_aucs.mean(), np.percentile(fixed_aucs,2.5), np.percentile(fixed_aucs,97.5)))
from collections import Counter
print('most frequent chosen weights:', Counter(chosen).most_common(3))
