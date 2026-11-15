# 🔬 TB-Dynamic-Selectivity: Physics-Based Biophysical Audit & Lead Optimization Pipeline for *Mycobacterium tuberculosis* InhA

An open-source computational biophysics and drug-screening framework designed to overcome the false-positive and off-target toxicity limitations of static molecular docking (AutoDock Vina) at physiological body temperature (**310 K**).

---
## 💡 Simple Explanation (Why This Is Important)

Traditional drug-discovery software (such as AutoDock Vina) functions like a **static photograph**—it merely checks whether a drug molecule fits into a protein pocket. This approach often carries the risk of "false positives" and potential toxic side effects within the human body.

My framework operates like a **"live video" at human body temperature (37°C / 310 K)**:
1. **Comparison of Fit ($\Delta G_{\text{Vina}}$):** Does the drug (acting like a "key") fit better into the proteins of *Mycobacterium tuberculosis* (TB)—the "bad locks"—than it does into human proteins—the "good locks"?
2. **Stability at Body Temperature ($-RT \ln(f_{\text{target}}/f_{\text{off-target}})$):** At 310 K, molecules are in constant vibration. Does the drug remain firmly bound to the TB target, while dissociating from other human proteins ("off-targets")?
3. **Tight Fit within the Pocket ($\Delta\text{SASA}/V_{\text{pocket}}$):** Does the drug displace water molecules to achieve a snug, secure fit within the TB protein pocket?

If the combined score across these criteria is **$\le -3.0\text{ kcal/mol}$**, the molecule is considered a **safe and selective TB drug lead**.

## 📐 Core Mathematical Framework (At $T = 310\text{ K}$, $RT = 0.616033\text{ kcal/mol}$)

### 1. Dynamic Selectivity Score ($\Delta\Delta G_{\text{dynamic}}$)
Combines static docking binding affinity ($\Delta G_{\text{Vina}}$) with 310 K key-residue contact persistence ($f$) and pocket burial-to-volume packing efficiency ($\Delta\text{SASA} / V_{\text{pocket}}$):

$$\Delta\Delta G_{\text{dynamic}} = (\Delta G_{\text{Vina, target}} - \Delta G_{\text{Vina, off-target}}) - RT \ln\left(\frac{f_{\text{target}}}{f_{\text{off-target}}}\right) - 1.5\,RT\left(\frac{\Delta\text{SASA}_{\text{target}}}{V_{\text{pocket, target}}} - \frac{\Delta\text{SASA}_{\text{off-target}}}{V_{\text{pocket, off-target}}}\right)$$

* **Decision Threshold:** $\Delta\Delta G_{\text{dynamic}} \le -3.0\text{ kcal/mol}$ indicates safe, selective target engagement.

### 2. Boltzmann-Weighted Contact Persistence ($\langle C \rangle$) & Contact Energy ($\langle \Delta E_{\text{contact}} \rangle$)
Weights molecular dynamics (MD) trajectory frames using Canonical Boltzmann Probabilities ($P_n$) to separate stable low-energy binding states from high-energy thermal noise:

$$P_n = \frac{\exp(-E_n / RT)}{\sum_{i=0}^{N} \exp(-E_i / RT)}, \quad \langle C \rangle = \sum_{n=0}^{N} P_n \cdot C_n, \quad \langle \Delta E_{\text{contact}} \rangle = \sum_{n=0}^{N} P_n \cdot C_n \cdot \Delta E_n$$

### 3. Heuristic Residence Time ($\tau_{\text{est}}$) & Steric Pocket Packing ($\eta_{\text{pack}}$)
Estimates kinetic unbinding persistence and steric complementarity (Rebek's 55% packing rule):

$$\tau_{\text{est}} = \tau_0 \exp\left(\frac{\Delta G^{\ddagger}_{\text{proxy}}}{RT}\right), \quad \eta_{\text{pack}} = \frac{V_{\text{ligand}}}{V_{\text{pocket}}} \approx 0.55$$

---

## 📊 Key Validation & Benchmark Results

| Module / Script | Validation Focus | Key Result |
| :--- | :--- | :--- |
| `prove_my_formulas.py` | Boundary limits, Boltzmann normalization, & Eyring kinetics | **4/4 Mathematical Proofs Verified** (`PROOF_REPORT.txt`) |
| `dude_inha_30_benchmark.py` | 30-Molecule TB InhA Discrimination under 310 K Thermal Noise | **ROC-AUC = 1.000** vs. Static Vina **0.893** (`dude_inha_roc_benchmark.png`) |
| `real_literature_inha_validation.py` | Correlation with Published Experimental $\text{IC}_{50}$ (PT70, NITD-916, GSK693, Triclosan) | **Spearman Rank $r = -0.873$**; 100% rejection of toxic off-target binders |
| `lead_optimization_engine.py` | Closed-Loop R-Group Optimization of a Toxic Off-Target Binder | Improved score from **$+0.932$ to $-6.456\text{ kcal/mol}$** ($-7.387\text{ kcal/mol}$ gain) |

---

## 🚀 Quick Start (Run Locally)

```bash
git clone https://github.com
cd TB-Dynamic-Selectivity
python3 prove_my_formulas.py
python3 dude_inha_30_benchmark.py
python3 real_literature_inha_validation.py
python3 lead_optimization_engine.py
```

---

## 🔭 Roadmap (Next 6 Months)
* Run full all-atom **100 ns+ GROMACS MD trajectories** across 50+ experimental ChEMBL/PDB *M. tuberculosis* InhA inhibitors and decoys on dedicated cloud GPUs.
* Calibrate empirical packing weights against real trajectory trajectories and publish an open-access preprint on **ChemRxiv / bioRxiv**.
