import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdDetermineBonds
from rdkit.Geometry import Point3D
TARGETS = {
 '4OIM': ('CCCCCCc1ccc(Oc2ccccc2C#N)c(O)c1', 'lig_4OIM_A.pdb', 'redock_4OIM.pdbqt'),
 '2X23': ('CCCCCCc1ccc(Oc2ccccc2C)c(O)c1', 'tcu_A.pdb', 'tcu_redock_out.pdbqt'),
 '4R9S': ('CC1(C)CCC(Cc2cc(O)c(-c3ccccc3)c(=O)[nH]2)CC1', 'lig_4R9S_A.pdb', 'redock_4R9S.pdbqt'),
 '5JFO': ('Cc1ccn([C@@H](C)c2nnc(Nc3ccn(Cc4c(F)cccc4Cl)n3)s2)n1', 'lig_5JFO_A.pdb', 'redock_5JFO.pdbqt'),
 '2NSD': ('Cc1ccc(cc1)C(=O)N2CCC(CC2)Cc3ccccc3', 'lig_2NSD_A.pdb', 'redock_2NSD.pdbqt'),
 '1P45': ('Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl', 'lig_1P45_A.pdb', 'redock_1P45.pdbqt'),
}
def elem(t):
    if t.upper()=='CL': return 'Cl'
    if t in ('A',): return 'C'
    if t=='OA': return 'O'
    if t=='NA': return 'N'
    return t[0].upper() if t[0].upper()+t[1:].lower() not in ('Cl',) else 'Cl'
def coords(fn, model=None):
    el=[];xyz=[];mdl=None
    for l in open(fn):
        if l.startswith('MODEL'): mdl=int(l.split()[1])
        if model is not None and mdl!=model: continue
        if model is None and not l.startswith('HETATM'): continue
        if not l.startswith(('ATOM','HETATM')): continue
        if fn.endswith('pdbqt'):
            t=l[77:79].strip()
            if t=='HD': continue
            e=elem(t)
        else:
            t=(l[76:78].strip() or l[12:14].strip())
            if t=='H': continue
            e=t.capitalize() if t.upper()=='CL' else t[0].upper()
        el.append(e); xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    m=Chem.RWMol()
    for e in el: m.AddAtom(Chem.Atom(e))
    conf=Chem.Conformer(len(el))
    for i,p in enumerate(xyz): conf.SetAtomPosition(i, Point3D(*p))
    m.AddConformer(conf, assignId=True)
    return m.GetMol()
for pid,(smi,cf,df) in TARGETS.items():
    try:
        tmpl=Chem.MolFromSmiles(smi)
        cm=coords(cf); dm=coords(df, model=1)
        rdDetermineBonds.DetermineConnectivity(cm); rdDetermineBonds.DetermineConnectivity(dm)
        cm2=AllChem.AssignBondOrdersFromTemplate(tmpl, cm)
        dm2=AllChem.AssignBondOrdersFromTemplate(tmpl, dm)
        print(f'{pid} | heavy crystal={cm.GetNumAtoms()} docked={dm.GetNumAtoms()} | CalcRMS = {rdMolAlign.CalcRMS(dm2, cm2):.2f} A', flush=True)
    except Exception as e:
        print(pid, 'ERROR', str(e)[:120], flush=True)
