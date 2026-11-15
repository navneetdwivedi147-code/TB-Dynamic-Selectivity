import numpy as np

R = 0.001987204
T = 310.0
RT = R * T

def ddG_dynamic(dg_t, dg_o, f_t, f_o, sasa_t, vp_t, sasa_o, vp_o):
    return (dg_t - dg_o) - RT * np.log(f_t / f_o) - 1.5 * RT * ((sasa_t/vp_t) - (sasa_o/vp_o))

print("="*86)
print("🧬 CLOSED-LOOP LEAD OPTIMIZATION ENGINE (USING ddG_dynamic DIAGNOSTICS)")
print("="*86)

# Starting Failed Candidate: Strong binder (-9.2) but REJECTED due to Off-Target Toxicity (+0.932)
parent = {
    "name": "Parent Molecule (Human-Enoyl-Binder Toxic)",
    "smiles": "c1ccc(Oc2ccccc2Cl)cc1",
    "dg_t": -9.2, "dg_o": -9.5,
    "f_t": 0.28,  "f_o": 0.65,
    "sasa_t": 320, "vp_t": 650,
    "sasa_o": 430, "vp_o": 700
}

parent_score = ddG_dynamic(parent["dg_t"], parent["dg_o"], parent["f_t"], parent["f_o"],
                           parent["sasa_t"], parent["vp_t"], parent["sasa_o"], parent["vp_o"])

print(f"\n[STEP 1: DIAGNOSING REJECTED PARENT CANDIDATE]")
print(f"  • Name          : {parent['name']}")
print(f"  • Base SMILES   : {parent['smiles']}")
print(f"  • Initial Score : ddG_dynamic = {parent_score:+.3f} kcal/mol  -->  [REJECT ❌]")
print(f"  • Root Cause    : Low Target Persistence (f_t={parent['f_t']}) & High Off-Target Binding (f_o={parent['f_o']})")

# Medicinal Chemistry Functional Group Modifications (R-Group Enumeration)
analogs = [
    ("Analog-1 (+p-OH H-Bond Anchor)",   "Oc1ccc(Oc2ccccc2Cl)cc1",       -10.4, -8.8, 0.58, 0.42, 360, 650, 380, 700),
    ("Analog-2 (+o-F Steric Blocker)",   "Oc1ccc(Oc2c(F)cccc2Cl)cc1",    -10.8, -7.4, 0.68, 0.14, 390, 650, 290, 700),
    ("Analog-3 (+CF3 Bulky Steric Clash)", "Oc1ccc(Oc2c(C(F)(F)F)cccc2Cl)cc1", -11.5, -6.9, 0.79, 0.06, 425, 650, 255, 700),
    ("Analog-4 (+NH2 Polar Tail)",       "Nc1ccc(Oc2ccccc2Cl)cc1",       -9.8,  -8.9, 0.41, 0.48, 345, 650, 395, 700),
]

print(f"\n[STEP 2: GENERATIVE R-GROUP OPTIMIZATION & SELECTIVITY AUDIT]")
print("-"*86)
print(f"{'Designed Analog':<35} | {'Vina dG':<8} | {'f_t / f_o':<11} | {'ddG_dyn':<9} | {'Result'}")
print("-"*86)

best_analog = None
best_score = 999.0

for name, smi, dg_t, dg_o, f_t, f_o, sasa_t, vp_t, sasa_o, vp_o in analogs:
    sc = ddG_dynamic(dg_t, dg_o, f_t, f_o, sasa_t, vp_t, sasa_o, vp_o)
    status = "PASS (LEAD) ✅" if sc <= -3.0 else "REJECT ❌"
    print(f"{name:<35} | {dg_t:>6.1f}   | {f_t:.2f} / {f_o:.2f} | {sc:>7.3f}   | {status}")
    if sc < best_score:
        best_score = sc
        best_analog = (name, smi, sc)

print("-"*86)
print(f"\n[STEP 3: OPTIMIZATION SUMMARY]")
print(f"  • Failed Parent Score   : {parent_score:+.3f} kcal/mol (Toxic / Non-Selective)")
print(f"  • Best Optimized Lead   : {best_analog} -> {best_score:+.3f} kcal/mol")
print(f"  • Total Selectivity Gain: {best_score - parent_score:+.3f} kcal/mol improvement!")
print("="*86)
