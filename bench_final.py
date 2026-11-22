import subprocess, numpy as np
from sklearn.metrics import roc_auc_score
from v4_ensemble import run, S_v4, E_CUT
import freesasa
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
    return [float(x) for x in subprocess.check_output("grep '^   [0-9]' " + logf + " | awk '{print $2}'", shell=True, text=True).split()]
E641 = [-9.883,-9.611,-9.025,-8.831,-8.654,-8.625,-8.042,-7.990,-7.973,-7.877]
dE_o, f_o, b_o, _, _ = run('641', 'lig_641_docked.pdbqt', 'inha_receptor.pdbqt', E641, '641')
rows=[]
for n,pout,plog in [('nitd916','nitd916_out.pdbqt','nitd916_vina.log'),('pt119','pt119_out.pdbqt','pt119_vina.log'),('gsk625','gsk625_out.pdbqt','gsk625_vina.log'),('ligand2nsd','ligand2nsd_out.pdbqt','ligand2nsd_vina.log'),('pt70','pt70_out.pdbqt','pt70_vina.log'),('triclosan','triclosan_out.pdbqt','triclosan_vina.log')]:
    E = vina_list(plog)
    dE,f,b,_,_ = run(n, pout, 'inha_receptor.pdbqt', E, n)
    rows.append((n,1,E[0],S_v4(dE,f,b,dE_o,f_o,b_o,free_sasa(pout))))
for line in open('new_actives_scores.txt'):
    p=line.split()
    if len(p)!=2: continue
    cid=p[0]
    E = vina_list('newact/%s_vina.log' % cid)
    dE,f,b,_,_ = run(cid, 'newact/%s_out.pdbqt' % cid, 'inha_receptor.pdbqt', E, cid)
    rows.append((cid,1,E[0],S_v4(dE,f,b,dE_o,f_o,b_o,free_sasa('newact/%s_out.pdbqt' % cid))))
for zid in [l.split()[1] for l in open('decoys_filtered.smi') if l.strip()]:
    E = vina_list('decoys/%s_vina.log' % zid)
    dE,f,b,_,_ = run(zid, 'decoys/%s_out.pdbqt' % zid, 'inha_receptor.pdbqt', E, zid)
    rows.append((zid,0,E[0],S_v4(dE,f,b,dE_o,f_o,b_o,free_sasa('decoys/%s_out.pdbqt' % zid))))
y=np.array([r[1] for r in rows]); v=np.array([r[2] for r in rows]); s=np.array([r[3] for r in rows])
gated=(s>=1e5).sum()
print('E_CUT=%.1f | actives=%d decoys=%d | gate-failed (S=1e6): %d' % (E_CUT, y.sum(), (y==0).sum(), gated))
print('AUC Vina = %.3f | AUC S_v4 = %.3f' % (roc_auc_score(y,-v), roc_auc_score(y,-s)))
rng=np.random.default_rng(42); bv=[]; bs=[]
for _ in range(2000):
    idx=rng.integers(0,len(y),len(y))
    if y[idx].sum()==0 or (1-y[idx]).sum()==0: continue
    bv.append(roc_auc_score(y[idx],-v[idx])); bs.append(roc_auc_score(y[idx],-s[idx]))
print('bootstrap 95%% CI  Vina [%.3f, %.3f]  S_v4 [%.3f, %.3f]' % (np.percentile(bv,2.5),np.percentile(bv,97.5),np.percentile(bs,2.5),np.percentile(bs,97.5)))
