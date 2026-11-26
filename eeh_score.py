# EEH (5G0S): stereo-correct build from RCSB InChI, docking log -> ensemble -> v4.0 S vs reference 641
import json, urllib.request
import numpy as np, freesasa
from rdkit import Chem
from rdkit.Chem import AllChem
from v4_ensemble import run, S_v4
freesasa.setVerbosity(freesasa.silent)
def build_sdf(out='structured33/EEH.sdf'):
    cc=json.load(urllib.request.urlopen('https://data.rcsb.org/rest/v1/core/chemcomp/EEH'))
    m=Chem.MolFromInchi(cc['rcsb_chem_comp_descriptor']['InChI'])
    mh=Chem.AddHs(m)
    AllChem.EmbedMolecule(mh, randomSeed=42); AllChem.MMFFOptimizeMolecule(mh)
    Chem.MolToMolFile(mh, out)
    return out
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
if __name__=='__main__':
    E641=[-9.883,-9.611,-9.025,-8.831,-8.654,-8.625,-8.042,-7.990,-7.973,-7.877]
    dE_o,f_o,b_o,_,_=run('641','lig_641_docked.pdbqt','inha_receptor.pdbqt',E641,'641')
    E=vina_list('structured33/EEH_vina.log')
    dE,f,b,_,_=run('EEH','structured33/EEH_out.pdbqt','inha_receptor.pdbqt',E,'EEH')
    sl=free_sasa('structured33/EEH_out.pdbqt')
    S=S_v4(dE,f,b,dE_o,f_o,b_o,sl)
    print('EEH | <dE>=%.3f <f_key>=%.4f <buried>=%.1f | freeSASA=%.1f | S=%.3f' % (dE,f,b,sl,S))
