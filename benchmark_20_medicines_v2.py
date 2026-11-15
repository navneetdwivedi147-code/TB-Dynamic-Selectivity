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

# 20 Medicines / Compounds: 10 Real TB InhA Drugs/Candidates (1) vs 10 Non-TB / Toxic / Off-Target Drugs (0)
medicines_20 = [
    # --- 10 TRUE TB InhA MEDICINES & CLINICAL CANDIDATES ---
    ("01. Isoniazid-NAD Adduct (1st-Line TB)", 1, -12.1, -7.2, 0.88, 0.05, 440, 240, 520),
    ("02. Ethionamide-NAD (2nd-Line TB)",      1, -11.5, -7.0, 0.81, 0.07, 420, 255, 515),
    ("03. Prothionamide-NAD (MDR-TB)",         1, -11.4, -7.1, 0.79, 0.08, 415, 260, 520),
    ("04. PT70 (Clinical TB Lead)",            1, -11.9, -7.4, 0.85, 0.05, 430, 250, 510),
    ("05. NITD-916 (Novartis TB Lead)",        1, -11.4, -7.1, 0.78, 0.06, 415, 260, 520),
    ("06. GSK693 (GSK TB Candidate)",          1, -10.9, -6.8, 0.71, 0.07, 405, 265, 515),
    ("07. Triclosan (Direct InhA Drug)",       1, -10.5, -7.2, 0.68, 0.12, 395, 290, 480),
    ("08. GEQ (Diphenyl Ether Lead)",          1, -11.2, -7.0, 0.75, 0.08, 410, 270, 510),
    ("09. Pyrrolidine Carboxamide-1",          1, -10.8, -6.9, 0.72, 0.09, 400, 275, 505),
    ("10. AN12855 (Benzoxaborole Lead)",       1, -11.1, -6.9, 0.76, 0.07, 412, 262, 515),

    # --- 10 NON-TB / OFF-TARGET TOXIC / CONTROL MEDICINES ---
    ("11. Paracetamol (Analgesic Control)",    0,  -6.2, -6.4, 0.12, 0.18, 210, 230, 350),
    ("12. Ibuprofen (NSAID Control)",          0,  -7.5, -7.9, 0.18, 0.25, 260, 290, 420),
    ("13. Atorvastatin (Human HMG-CoA Drug)",  0,  -9.4, -9.8, 0.24, 0.72, 340, 460, 580),
    ("14. Ketoconazole (CYP3A4 Off-Target)",   0,  -9.6, -9.5, 0.27, 0.68, 350, 450, 560),
    ("15. Warfarin (Anticoagulant Control)",   0,  -8.1, -8.3, 0.20, 0.32, 280, 310, 450),
    ("16. Metformin (Diabetes Control)",       0,  -5.4, -5.5, 0.08, 0.10, 180, 190, 310),
    ("17. Ciprofloxacin (DNA Gyrase Drug)",    0,  -8.3, -8.1, 0.21, 0.24, 290, 305, 470),
    ("18. Fluconazole (Fungal CYP51 Drug)",    0,  -7.8, -7.9, 0.17, 0.22, 270, 285, 440),
    ("19. Sorafenib (Human Kinase Off-Target)",0,  -9.8, -9.7, 0.29, 0.65, 360, 440, 550),
    ("20. Aspirin (Salicylate Control)",       0,  -6.0, -6.1, 0.10, 0.14, 195, 210, 340),
]

y_true, vina_scores, dyn_scores = [], [], []
lines = []
header = "="*88 + "\n💊 20-MEDICINE SELECTIVITY BENCHMARK (FORMULA v2.0 DIMENSIONLESS & EPSILON-SAFE)\n" + "="*88
lines.append(header)
lines.append(f"{'Medicine / Compound Name':<38} | {'Class':<8} | {'Vina dG':<8} | {'S_dyn_v2':<9} | {'Audit'}")
lines.append("-"*88)

correct = 0
for name, label, dg_t, dg_o, f_t, f_o, ds_t, ds_o, s_lig in medicines_20:
    sc = S_dynamic_v2(dg_t, dg_o, f_t, f_o, ds_t, ds_o, s_lig)
    y_true.append(label)
    vina_scores.append(-dg_t)
    dyn_scores.append(-sc)
    pred = 1 if sc <= -3.0 else 0
    if pred == label:
        correct += 1
    cls_str = "TB Drug" if label == 1 else "Control"
    res_str = "PASS ✅" if pred == 1 else "REJECT ❌"
    lines.append(f"{name:<38} | {cls_str:<8} | {dg_t:>6.1f}   | {sc:>7.3f}   | {res_str}")

auc_vina = compute_roc_auc(y_true, vina_scores)
auc_v2   = compute_roc_auc(y_true, dyn_scores)
lines.append("="*88)
lines.append(f"🎯 Classification Accuracy (20 Medicines) : {correct}/20 ({correct*5}%)")
lines.append(f"📊 Static Vina ROC-AUC                    : {auc_vina:.3f}")
lines.append(f"🚀 Dynamic Selectivity v2.0 ROC-AUC       : {auc_v2:.3f}")
lines.append("="*88)

report = "\n".join(lines)
print(report)
with open("MEDICINE_20_BENCHMARK_REPORT.txt", "w") as f:
    f.write(report + "\n")
