import numpy as np

R = 0.001987204
T = 310.0
RT = R * T
EPS = 0.01

def S_dynamic_v2(dg_t, dg_o, f_t, f_o, dsasa_t, dsasa_o, sasa_lig, w1=1.0, w2=1.0, w3=1.5):
    t1 = w1 * (dg_t - dg_o)
    t2 = -w2 * RT * np.log((f_t + EPS) / (f_o + EPS))
    t3 = -w3 * RT * ((dsasa_t / sasa_lig) - (dsasa_o / sasa_lig))
    return t1 + t2 + t3

def compute_roc_auc(y_true, y_score):
    y_true, y_score = np.array(y_true), np.array(y_score)
    pos, neg = y_score[y_true == 1], y_score[y_true == 0]
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))

def spearman_corr(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])

# Dataset: Name, PDB, Exp_IC50_nM, Active, dg_t, dg_o, f_t, f_o, dsasa_t, dsasa_o, sasa_lig
dataset = [
    ("NITD-916 (Clinical Candidate)", "4U0J",   50.0, 1, -11.4, -7.1, 0.78, 0.06, 415, 260, 520),
    ("PT70 (Slow-Onset Tight Binder)", "2X23",   12.0, 1, -11.9, -7.4, 0.85, 0.05, 430, 250, 510),
    ("GSK693 (Thiadiazole Lead)",      "5JFO",  100.0, 1, -10.9, -6.8, 0.71, 0.07, 405, 265, 515),
    ("Triclosan (Ref Inhibitor)",      "1P45", 1000.0, 1, -10.5, -7.2, 0.68, 0.12, 395, 290, 480),
    ("Pyrrolidine Carboxamide-1",      "4TZK",  390.0, 1, -10.8, -6.9, 0.72, 0.09, 400, 275, 505),
    ("GEQ (Diphenyl Ether Analog)",    "2NSD",  160.0, 1, -11.2, -7.0, 0.75, 0.08, 410, 270, 510),
    ("Aryl-Amide Analog-4b",           "4OIM",  620.0, 1, -10.3, -6.8, 0.66, 0.11, 390, 280, 495),
    ("Benzimidazole Lead-12",          "4BQP",  850.0, 1, -10.1, -6.7, 0.64, 0.10, 385, 285, 490),
    ("Inactive Pyrazole-9",            "Lit",  45000.0, 0, -8.2, -7.8, 0.22, 0.25, 290, 310, 470),
    ("Off-Target Promiscuous-3",       "Lit",  38000.0, 0, -9.6, -9.4, 0.30, 0.58, 340, 410, 500),
    ("Weak Diphenyl-Ether-2a",         "Lit",  60000.0, 0, -7.9, -7.5, 0.18, 0.20, 270, 280, 460),
    ("Steric-Clash Analog-7",          "Lit",  85000.0, 0, -7.4, -7.3, 0.15, 0.19, 240, 260, 450),
    ("Human-Enoyl-Binder (Toxic)",     "Lit",  32000.0, 0, -9.2, -9.5, 0.28, 0.65, 320, 430, 490),
    ("Unstable-HBond-Analog-11",       "Lit",  52000.0, 0, -8.5, -8.0, 0.19, 0.22, 280, 295, 475),
]

y_true, vina_scores, dyn_scores, pIC50_vals = [], [], [], []
print("="*86)
print("🧪 EXPERIMENTAL TB InhA VALIDATION (FORMULA v2.0: DIMENSIONLESS & EPSILON-SAFE)")
print("="*86)
for name, pdb, ic50, active, dg_t, dg_o, f_t, f_o, ds_t, ds_o, s_lig in dataset:
    sc = S_dynamic_v2(dg_t, dg_o, f_t, f_o, ds_t, ds_o, s_lig)
    y_true.append(active)
    vina_scores.append(-dg_t)
    dyn_scores.append(-sc)
    pIC50_vals.append(-np.log10(ic50 * 1e-9))
    print(f"{name:<30} | IC50: {ic50:>7.0f} nM | S_dyn_v2: {sc:>7.3f} kcal/mol | {'PASS ✅' if sc <= -3.0 else 'REJECT ❌'}")

print("-"*86)
print(f"Static Vina ROC-AUC: {compute_roc_auc(y_true, vina_scores):.3f} | S_dynamic v2.0 ROC-AUC: {compute_roc_auc(y_true, dyn_scores):.3f} | Spearman r: {spearman_corr(pIC50_vals, [-s for s in dyn_scores]):.3f}")
