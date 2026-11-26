import numpy as np, freesasa
from sklearn.metrics import roc_auc_score
from rdkit import Chem
from rdkit.Chem import Descriptors
import v4_ensemble as V
freesasa.setVerbosity(freesasa.silent)
def elem(t):
    t=t.strip()
    if t in ('A','C'): return 'C'
    if t.startswith('OA') or t=='O': return 'O'
    if t.startswith('NA') or t=='N': return 'N'
    if t.startswith('SA') or t=='S': return 'S'
    if t in ('F','Cl','CL'): return 'CL' if t!='F' else 'F'
    if t.startswith('H'): return 'H'
    return t[0]
def free_sasa(fn):
    s=freesasa.Structure(); i=0; m=None
    for l in open(fn):
        if l.startswith('MODEL'): m=int(l.split()[1])
        if m==1 and l.startswith(('ATOM','HETATM')):
            t=l[77:79]; e = elem(t) if t.strip() else elem(l[12:14])
            if e!='H':
                i+=1; s.addAtom(e,'U',str(i),'A',float(l[30:38]),float(l[38:46]),float(l[46:54]))
    return freesasa.calc(s).totalArea()
def vina_list(logf):
    out=[]
    for l in open(logf, errors='ignore'):
        p=l.split()
        if len(p)>=2 and p[0].isdigit() and l.lstrip().startswith(tuple(str(k) for k in range(1,10))):
            try: out.append(float(p[1]))
            except: pass
    return out
def S_eps(dE,f,b,dEo,fo,bo,sl,eps):
    if dE>0.0: return 1e6
    RT=V.RT
    return (V.W1*(dE-dEo) - V.W2*RT*np.log((f+eps)/(fo+eps)) - V.W3*RT*(b/sl - bo/sl))
E641=[-9.883,-9.611,-9.025,-8.831,-8.654,-8.625,-8.042,-7.990,-7.973,-7.877]
dEo,fo,bo,_,_=V.run('641','lig_641_docked.pdbqt','inha_receptor.pdbqt',E641,'641')
slo=None
Y=[];rows=[]
for l in open('structured33_v4.txt'):
    p=l.split()
    if len(p)<8 or p[0]=='641' or 'NO_ENSEMBLE' in l: continue
    name=p[0]
    E=vina_list(f'structured33/{name}_vina.log')
    dE,f,b,_,_=V.run(name,f'structured33/{name}_out.pdbqt','inha_receptor.pdbqt',E,name)
    sl=free_sasa(f'structured33/{name}_out.pdbqt')
    rows.append(('act',name,E[0],dE,f,b,sl)); Y.append(1)
for l in list(open('decoy_scores_by_zinc.tsv'))[1:]:
    p=l.rstrip('\n').split('\t')
    if len(p)<4 or p[2]=='NA': continue
    zid=p[0]
    E=vina_list(f'decoys/{zid}_vina.log')
    if len(E)<8: continue
    dE,f,b,_,_=V.run(zid,f'decoys/{zid}_out.pdbqt','inha_receptor.pdbqt',E,zid)
    sl=free_sasa(f'decoys/{zid}_out.pdbqt')
    rows.append(('dec',zid,E[0],dE,f,b,sl)); Y.append(0)
print('actives', sum(Y), '| decoys', len(Y)-sum(Y))
np.save('eps_inputs.npy', np.array([[r[2],r[3],r[4],r[5],r[6]] for r in rows]))
np.save('eps_labels.npy', np.array(Y))
y=np.array(Y)
v=np.array([r[2] for r in rows])
for eps in (0.001,0.02,0.05):
    s=np.array([S_eps(r[3],r[4],r[5],dEo,fo,bo,r[6],eps) for r in rows])
    print('eps=%.3f | AUC S = %.3f | AUC Vina = %.3f' % (eps, roc_auc_score(y,-s), roc_auc_score(y,-v)))
