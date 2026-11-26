# SUPERSEDED: graph-only atom mapping gave wrong RMSD for ring/symmetric atoms. Use redock_rmsd_calc.py (RDKit CalcRMS, bond-aware).
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D
TARGETS = {
 '4R9S': ('CC1(C)CCC(Cc2cc(O)c(-c3ccccc3)c(=O)[nH]2)CC1', 'lig_4R9S_A.pdb', 'redock_4R9S.pdbqt'),
 '5JFO': ('Cc1ccn([C@@H](C)c2nnc(Nc3ccn(Cc4c(F)cccc4Cl)n3)s2)n1', 'lig_5JFO_A.pdb', 'redock_5JFO.pdbqt'),
 '2NSD': ('Cc1ccc(cc1)C(=O)N2CCC(CC2)Cc3ccccc3', 'lig_2NSD_A.pdb', 'redock_2NSD.pdbqt'),
 '1P45': ('Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl', 'lig_1P45_A.pdb', 'redock_1P45.pdbqt'),
 '4OIM': ('CCCCCCc1ccc(Oc2ccccc2C#N)c(O)c1', 'lig_4OIM_A.pdb', 'redock_4OIM.pdbqt'),
 '2X23': ('CCCCCCc1ccc(Oc2ccccc2C)c(O)c1', 'lig_2X23_A.pdb', 'redock_2X23.pdbqt'),
}
def graph_only(m):
    g=Chem.RWMol(m)
    for b in g.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    for a in g.GetAtoms(): a.SetIsAromatic(False); a.SetNoImplicit(True)
    return g.GetMol()
def build(el, xyz):
    m=Chem.RWMol()
    for e in el: m.AddAtom(Chem.Atom(e))
    conf=Chem.Conformer(len(el))
    for i,p in enumerate(xyz): conf.SetAtomPosition(i, Point3D(*p))
    m.AddConformer(conf, assignId=True)
    mm=m.GetMol(); rdDetermineBonds.DetermineConnectivity(mm); return mm
def crystal(fn):
    el=[]; xyz=[]
    for l in open(fn):
        if not l.startswith('HETATM'): continue
        e=(l[76:78].strip() or l[12:14].strip()).capitalize()
        if e=='H': continue
        el.append(e); xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return el, np.array(xyz)
def docked(fn):
    el=[]; xyz=[]; mdl=None
    for l in open(fn):
        if l.startswith('MODEL'): mdl=int(l.split()[1])
        if mdl!=1 or not l.startswith(('ATOM','HETATM')): continue
        t=l[77:79].strip()
        if t=='HD': continue
        el.append({'A':'C','OA':'O','NA':'N','Cl':'Cl','C':'C','N':'N','O':'O','S':'S'}.get(t,t[0]))
        xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return el, np.array(xyz)
def rmsd_sym(smi, cfile, dfile):
    rg=graph_only(Chem.MolFromSmiles(smi))
    cel,cx=crystal(cfile); del_,dx=docked(dfile)
    cg=graph_only(build(cel,cx)); dg=graph_only(build(del_,dx))
    mc=rg.GetSubstructMatch(cg)
    if len(mc)!=len(cel): return None, 'crystal map fail'
    mds=dg.GetSubstructMatches(rg, uniquify=False, maxMatches=5000)
    best=None
    for md in mds:
        if len(md)!=len(del_): continue
        cp=np.zeros((len(cel),3)); dp=np.zeros((len(cel),3))
        for i,si in enumerate(mc): cp[si]=cx[i]
        for i,si in enumerate(md): dp[si]=dx[i]
        r=np.sqrt(((cp-dp)**2).sum(1).mean())
        if best is None or r<best: best=r
    return best, '%d heavy, %d matches' % (len(cel), len(mds))
if __name__=='__main__':
    for pid,(smi,cfile,dfile) in TARGETS.items():
        try:
            r,info=rmsd_sym(smi,cfile,dfile)
            print(pid, '| RMSD = %s | %s' % (('%.2f A' % r) if r is not None else 'NA', info))
        except Exception as e:
            print(pid, 'ERROR', e)
