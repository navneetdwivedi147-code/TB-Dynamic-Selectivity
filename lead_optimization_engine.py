import numpy as np

R = 0.001987204
T = 310.0
RT = R * T
EPS = 0.01

def S_dynamic_v2(dg_t, dg_o, f_t, f_o, dsasa_t, dsasa_o, sasa_lig, w1=1.0, w2=1.0, w3=1.5):
    return w1*(dg_t - dg_o) - w2*RT*np.log((f_t + EPS)/(f_o + EPS)) - w3*RT*((dsasa_t/sasa_lig) - (dsasa_o/sasa_lig))

parent_sc = S_dynamic_v2(-9.2, -9.5, 0.28, 0.65, 320, 430, 490)
analogs = [
    ("Analog-1 (+p-OH H-Bond Anchor)",     -10.4, -8.8, 0.58, 0.42, 360, 380, 505),
    ("Analog-2 (+o-F Steric Blocker)",     -10.8, -7.4, 0.68, 0.14, 390, 290, 515),
    ("Analog-3 (+CF3 Bulky Steric Clash)", -11.5, -6.9, 0.79, 0.06, 425, 255, 540),
    ("Analog-4 (+NH2 Polar Tail)",         -9.8,  -8.9, 0.41, 0.48, 345, 395, 505),
]

print("="*86)
print("🧬 CLOSED-LOOP LEAD OPTIMIZATION ENGINE (FORMULA v2.0 DIMENSIONLESS)")
print("="*86)
print(f"Parent Toxic Binder Score: {parent_sc:+.3f} kcal/mol [REJECT ❌]")
for name, dg_t, dg_o, f_t, f_o, ds_t, ds_o, s_lig in analogs:
    sc = S_dynamic_v2(dg_t, dg_o, f_t, f_o, ds_t, ds_o, s_lig)
    print(f"  • {name:<35} -> S_dyn_v2 = {sc:+.3f} kcal/mol [{'PASS (LEAD) ✅' if sc <= -3.0 else 'REJECT ❌'}]")
print("="*86)
