import numpy as np, freesasa
freesasa.setVerbosity(freesasa.silent)
RT = 0.001987204 * 310.0
W1, W2, W3, EPS, E_CUT, TAU0 = 1.0, 1.0, 1.5, 0.02, 0.0, 1.0
R = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'F':1.47,'CL':1.75}
def elem(t):
    t=t.strip()
    if t in ('A','C'): return 'C'
    if t.startswith('OA') or t=='O': return 'O'
    if t.startswith('NA') or t=='N': return 'N'
    if t.startswith('SA') or t=='S': return 'S'
    if t in ('F','Cl','CL'): return 'CL' if t!='F' else 'F'
    if t.startswith('H'): return 'H'
    return t[0]
def is_heavy(l):
    t=l[77:79]; return (elem(t) if t.strip() else elem(l[12:14])) != 'H'
def read_modes(fn):
    modes, cur, m = {}, [], None
    for l in open(fn):
        if l.startswith('MODEL'): m=int(l.split()[1]); cur=[]
        elif l.startswith('ENDMDL'): modes[m]=cur
        elif l.startswith(('ATOM','HETATM')) and m is not None and is_heavy(l): cur.append(l)
    return modes
def read_rec(fn):
    return [l for l in open(fn) if l.startswith(('ATOM','HETATM')) and is_heavy(l)]
def xyz(l): return np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
def el(l):
    t=l[77:79]; return elem(t) if t.strip() else elem(l[12:14])
def rad(l): return R.get(el(l),1.7)
def boltz(E):
    E=np.array(E); w=np.exp(-(E-E.min())/RT); return w/w.sum()
def pose_metrics(lig, rec):
    lx=np.array([xyz(l) for l in lig]); rx=np.array([xyz(l) for l in rec])
    tyr = [l for l in rec if l[17:20].strip()=='TYR' and l[22:26].strip()=='158']
    tx = np.array([xyz(l) for l in tyr])
    d=np.linalg.norm(lx[:,None,:]-tx[None,:,:],axis=2)
    dmin = d.min(axis=1)
    w = np.clip((4.2 - dmin) / (4.2 - 3.0), 0.0, 1.0)
    f_key = float(w.mean())
    near=[l for l in rec if np.min(np.linalg.norm(lx-xyz(l),axis=1))<10.0]
    s0=freesasa.Structure()
    for i,l in enumerate(lig):
        x,y,z=xyz(l); s0.addAtom(el(l),'U',str(i+1),'A',float(x),float(y),float(z))
    free=freesasa.calc(s0).totalArea()
    s=freesasa.Structure()
    for i,l in enumerate(lig+near):
        x,y,z=xyz(l); s.addAtom(el(l),'U',str(i+1),'A',float(x),float(y),float(z))
    bound=sum(freesasa.calc(s).atomArea(i) for i in range(len(lig)))
    return f_key, free-bound
def run(name, dock_file, rec_file, energies, pose_source):
    modes = read_modes(dock_file)
    rec = read_rec(rec_file)
    keys = sorted(modes)
    E = np.array(energies[:len(keys)])
    P = boltz(E)
    fk, bur = [], []
    for k in keys:
        f,b = pose_metrics(modes[k], rec)
        fk.append(f); bur.append(b)
    fk=np.array(fk); bur=np.array(bur)
    dE_ens = float(P@E)
    fk_ens = float(P@fk)
    bur_ens = float(P@bur)
    return dE_ens, fk_ens, bur_ens, E, P
def S_v4(dE_t, f_t, bur_t, dE_o, f_o, bur_o, sasa_lig):
    if dE_t > E_CUT: return 1e6
    return (W1*(dE_t-dE_o) - W2*RT*np.log((f_t+EPS)/(f_o+EPS)) - W3*RT*(bur_t/sasa_lig - bur_o/sasa_lig))
def tau_v4(S): return TAU0*np.exp(-S/RT)
if __name__ == '__main__':
    E_mtb = [-8.560,-8.493,-8.219,-8.161,-7.880,-7.563,-7.535,-7.274,-7.171]
    dE, fk, bur, E, P = run('MTB-NL4', 'mtb_nl4_out.pdbqt', 'inha_receptor.pdbqt', E_mtb, 'mtb')
    print('MTB-NL4 ensemble: <dE_int>=%.3f  <f_key>=%.4f  <buried>=%.1f A^2' % (dE, fk, bur))
    print('MTB-NL4 P_n:', np.round(P,3))
