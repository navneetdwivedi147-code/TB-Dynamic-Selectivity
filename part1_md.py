import os, subprocess, glob
os.chdir(os.path.expanduser("~/TB_MD_Simulation"))

for fname in ["InhA_receptor.pdb", "Novel_Lead_4_5CF3_pose1.mol2"]:
    if not os.path.exists(fname):
        found = glob.glob(f"/mnt/c/Users/*/Downloads/{fname}")
        if found: subprocess.run(f"cp '{found}' .", shell=True)
        else: raise FileNotFoundError(f"❌ {fname} nahi mili! Pehle Google Drive se {fname} download karein.")

def run(cmd):
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0: print(f"⚠️ Err on: {cmd[:60]}...\n{res.stderr[-250:]}")
    return res

print("⏳ 1/4 Preparing Receptor & GAFF2 Ligand Topology...")
run("rm -f topol.top posre.itp protein.gro complex.gro box.gro solv.gro solv_ions.gro em.* npt_50ps.* md_250ps.* '#*'")
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fixer = PDBFixer(filename="InhA_receptor.pdb")
fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.4)
with open("protein_clean.pdb", "w") as f: PDBFile.writeFile(fixer.topology, fixer.positions, f)

run("gmx pdb2gmx -f protein_clean.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh")
run("acpype -i Novel_Lead_4_5CF3_pose1.mol2 -b LIG -c bcc -n 0 -a gaff2")

with open("protein.gro") as f: p_lines = f.readlines()
with open("LIG.acpype/LIG_GMX.gro") as f: l_lines = f.readlines()
n_prot, n_lig = int(p_lines.strip()), int(l_lines.strip())
with open("complex.gro", "w") as f:
    f.write("InhA_MTB-NL4_Complex\n" + f"{n_prot + n_lig}\n")
    f.writelines(p_lines[2:-1]); f.writelines(l_lines[2:-1]); f.write(p_lines[-1])

with open("topol.top") as f: top = f.read()
top = top.replace('#include "amber99sb-ildn.ff/forcefield.itp"', '#include "amber99sb-ildn.ff/forcefield.itp"\n#include "LIG.acpype/LIG_GMX.itp"') + "\nLIG                 1\n"
with open("topol.top", "w") as f: f.write(top)

print("⏳ 2/4 Solvation, Ions & Energy Minimization...")
run("gmx editconf -f complex.gro -o box.gro -c -d 1.0 -bt dodecahedron")
run("gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top")
with open("ions.mdp", "w") as f: f.write("integrator = steep\nnsteps = 500\nemtol = 1000.0\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n")
run("gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2")
run("echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15")

with open("em.mdp", "w") as f: f.write("integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 5000\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n")
run("gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2")
run("gmx mdrun -v -deffnm em")

print("⏳ 3/4 Running 50 ps NPT Equilibration...")
with open("npt.mdp", "w") as f:
    f.write("integrator = md\ndt = 0.002\nnsteps = 25000\nnstxout-compressed = 5000\nnstenergy = 1000\ncontinuation = no\nconstraint_algorithm = lincs\nconstraints = h-bonds\ncutoff-scheme = Verlet\nrcoulomb = 1.0\nrvdw = 1.0\ncoulombtype = PME\ntcoupl = V-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 310\npcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\ncompressibility = 4.5e-5\npbc = xyz\ngen_vel = yes\ngen_temp = 310\n")
run("gmx grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt_50ps.tpr -maxwarn 2")
run("gmx mdrun -v -deffnm npt_50ps")

print("⏳ 4/4 Running 250 ps Pilot MD...")
with open("md_250ps.mdp", "w") as f:
    f.write("integrator = md\ndt = 0.002\nnsteps = 125000\nnstxout-compressed = 500\nnstenergy = 500\ncontinuation = yes\nconstraint_algorithm = lincs\nconstraints = h-bonds\ncutoff-scheme = Verlet\nrcoulomb = 1.0\nrvdw = 1.0\ncoulombtype = PME\ntcoupl = V-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 310\npcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\ncompressibility = 4.5e-5\npbc = xyz\ngen_vel = no\n")
run("gmx grompp -f md_250ps.mdp -c npt_50ps.gro -t npt_50ps.cpt -p topol.top -o md_250ps.tpr -maxwarn 2")
run("gmx mdrun -v -deffnm md_250ps")
print("✅ PART 1 COMPLETE! md_250ps.xtc exists:", os.path.exists("md_250ps.xtc"))
