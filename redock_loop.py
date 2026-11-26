import numpy as np, subprocess, os, sys
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D
TARGETS = [('4R9S','3KY','CC1(C)CCC(Cc2cc(O)c(-c3ccccc3)c(=O)[nH]2)CC1'),
           ('5JFO','6KA','Cc1ccn([C@@H](C)c2nnc(Nc3ccn(Cc4c(F)cccc4Cl)n3)s2)n1'),
           ('4OIM','JUS','CCCCCCc1ccc(Oc2ccccc2C#N)c(O)c1'),
           ('1P45','TCL','Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl'),
           ('2NSD','4PI','Cc1ccc(cc1)C(=O)N2CCC(CC2)Cc3ccccc3'),
           ('1P44','GEQ','O=C(N1CCN(CC1)C2c3ccccc3c4ccccc24)c5ccc6[nH]ccc6c5')]
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
for pid,h,smi in TARGETS:
    pdb=f'{pid}.pdb'; ch='A'
    prot=[l for l in open(pdb) if l.startswith('ATOM') and l[21]==ch]
    lig=[l for l in open(pdb) if l.startswith('HETATM') and l[17:20]==h and l[21]==ch]
    if not lig or not prot:
        print(pid,'SKIP no ligand/receptor chain A'); continue
    lxyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lig])
    c=lxyz.mean(0); ext=lxyz.max(0)-lxyz.min(0)
    size=np.round(ext+8,1)
    rec=f'rec_{pid}_A.pdb'; open(rec,'w').writelines(prot)
    subprocess.run(['obabel',rec,'-h','-xr','-O',f'rec_{pid}_A.pdbqt'],capture_output=True)
    lp=f'lig_{pid}_A.pdb'; open(lp,'w').writelines(lig)
    subprocess.run(['obabel',lp,'-O',f'lig_{pid}_A.sdf'],capture_output=True)
    subprocess.run(['obabel',f'lig_{pid}_A.sdf','-O',f'lig_{pid}_A.pdbqt','-h'],capture_output=True)
    out=f'redock_{pid}.pdbqt'; logf=f'redock_{pid}.log'
    subprocess.run(['vina','--receptor',f'rec_{pid}_A.pdbqt','--ligand',f'lig_{pid}_A.pdbqt',
        '--center_x',str(round(c[0],2)),'--center_y',str(round(c[1],2)),'--center_z',str(round(c[2],2)),
        '--size_x',str(size[0]),'--size_y',str(size[1]),'--size_z',str(size[2]),
        '--exhaustiveness','32','--num_modes','9','--out',out],stdout=open(logf,'w'),stderr=subprocess.STDOUT)
    score=None
    for l in open(logf, errors='ignore'):
        if l.lstrip().startswith('1 ') or l.startswith('   1 '):
            p=l.split()
            if len(p)>=2 and p[0]=='1': score=p[1]; break
    ref=Chem.MolFromSmiles(smi); rg=graph_only(ref)
    cel=[];cx=[]
    for l in lig:
        e=(l[76:78].strip() or l[12:14].strip()).capitalize()
        if e=='H': continue
        cel.append(e); cx.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    cx=np.array(cx)
    dels=[];dx=[];mdl=None
    if os.path.exists(out):
        for l in open(out):
            if l.startswith('MODEL'): mdl=int(l.split()[1])
            if mdl!=1 or not l.startswith(('ATOM','HETATM')): continue
            t=l[77:79].strip()
            if t=='HD': continue
            el={'A':'C','OA':'O','NA':'N','Cl':'Cl','C':'C','N':'N','O':'O','S':'S'}.get(t,t[0])
            dels.append(el); dx.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    dx=np.array(dx) if dx else np.zeros((0,3))
    if not dels:
        print(pid,'NO DOCKED POSE'); continue
    cg=graph_only(build(cel,cx)); dg=graph_only(build(dels,dx))
    mc=rg.GetSubstructMatch(cg); md=rg.GetSubstructMatch(dg)
    if len(mc)!=len(cel) or len(md)!=len(dels):
        print(pid,'score',score,'| MAP FAIL', len(mc),len(cel),len(md),len(dels)); continue
    cpos=np.zeros((len(cel),3)); dpos=np.zeros((len(dels),3))
    for i,si in enumerate(mc): cpos[si]=cx[i]
    for i,si in enumerate(md): dpos[si]=dx[i]
    rmsd=np.sqrt(((cpos-dpos)**2).sum(1).mean())
    print(f'{pid} ({h}) | vina1={score} | heavy={len(cel)} | size={tuple(size)} | RMSD={rmsd:.2f} A', flush=True)
