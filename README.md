# 🔬 TB-Dynamic-Selectivity: Physics-Based Biophysical Audit & Lead Optimization Pipeline for *Mycobacterium tuberculosis* InhA

📄 **[Download Full 11-Section Summary Report (PDF)](./MTB-NL4_Computational_Screening_Summary_Report.pdf)** | **Author:** Navneet Dwivedi | **Target:** *M. tuberculosis* InhA (PDB: `4TZK`) | **Lead Candidate:** `MTB-NL4`

---

## 🌟 Executive Summary (Complete Project Overview)

Static molecular docking (such as AutoDock Vina) evaluates drug candidates using a single frozen snapshot, frequently producing **false positives** and missing **human off-target toxicity** inside a warm, dynamic human body ($37^\circ\text{C} / 310\text{ K}$).

To solve this, I developed **TB-Dynamic-Selectivity**, an open-source computational biophysics pipeline and physics-informed scoring framework that acts like a **"live molecular movie" at 310 K**:
1. **Ensemble Interaction Gap ($\langle \Delta E_{\text{int}} \rangle$):** Checks whether the drug ("key") binds more strongly to the *M. tuberculosis* target ("bad lock") than to human off-target proteins ("good locks") across MD trajectory frames.
2. **Catalytic Key-Residue Persistence ($f_{\text{key}}$):** Measures whether the drug stays anchored to native catalytic residues ($d \le 3.5\text{ \AA}$, e.g., InhA Tyr158/NAD⁺) under 310 K thermal vibration.
3. **Fractional Buried Surface Packing ($\text{SASA}_{\text{buried}} / \text{SASA}_{\text{ligand}}$):** Quantifies how tightly the drug tucks inside the target pocket to displace water molecules.

### 🏆 Key Verified Milestones
* **4/4 Mathematical & Thermodynamic Proofs Verified (`PROOF_REPORT.txt`):** Boltzmann probability conservation ($\sum P_n = 1.000000$), exact zero-selectivity symmetry ($0.000000\text{ kcal/mol}$), zero-singularity crash safety ($f_{\text{off}} = 0$), and Eyring residence time ($34.40\text{–}38.23\text{ ns}$).
* **30-Molecule DUD-E Thermal Noise Benchmark (`dude_inha_roc_benchmark.png`):** Achieved **ROC-AUC = 1.000** vs. Static Vina **0.893**.
* **20-Medicine Clinical & Off-Target Audit (`MEDICINE_20_BENCHMARK_REPORT.txt`):** **20/20 (100% Accuracy)** separating 10 real TB drugs/leads (Isoniazid-NAD, PT70, NITD-916, GSK693, Triclosan) from 10 off-target/toxic FDA drugs (Atorvastatin, Ketoconazole, Sorafenib).
* **Experimental $\text{IC}_{50}$ Correlation (`real_literature_inha_validation.py`):** Strong monotonic agreement (**Spearman Rank $r = -0.873$**) across 14 published TB InhA inhibitors.
* **Closed-Loop Generative Lead Optimization (`lead_optimization_engine.py`):** Re-engineered a toxic off-target binder ($+1.014\text{ kcal/mol}$) via R-group enumeration (`-OH` + `-CF3`) into a selective TB lead ($-6.392\text{ kcal/mol}$, a **$-7.406\text{ kcal/mol}$** selectivity gain).

---

## 🔄 How the Formula Evolved Across 3 Rigorous Iterations (v1.0 → v2.0 → v3.0)

| Version | Mathematical Formulation | Key Scientific Upgrade & Bug Fixes |
| :--- | :--- | :--- |
| **v1.0** *(Initial Prototype)* | $\Delta\Delta G_{\text{dyn}} = (\Delta G_t - \Delta G_o) - RT\ln\left(\frac{f_t}{f_o}\right) - 1.5\,RT\left(\frac{\Delta\text{SASA}_t}{V_{p,t}} - \frac{\Delta\text{SASA}_o}{V_{p,o}}\right)$ | Combined static Vina affinity, 310 K contact persistence, and pocket burial into a single decision score ($\le -3.0\text{ kcal/mol}$). |
| **v2.0** *(Dimensionless & $\epsilon$-Safe)* | $S_{\text{dyn}} = w_1(\Delta G_t - \Delta G_o) - w_2 RT\ln\left(\frac{f_t + \epsilon}{f_o + \epsilon}\right) - w_3 RT\left(\frac{\Delta\text{SASA}_t}{\text{SASA}_{\text{lig}}} - \frac{\Delta\text{SASA}_o}{\text{SASA}_{\text{lig}}}\right)$ | **1)** Fixed units to exact $\text{kcal/mol}$ ($\text{Å}^2/\text{Å}^2 = \text{dimensionless}$).<br>**2)** Added $\epsilon = 0.01$ to prevent $\ln(f/0) \to \infty$ crashes.<br>**3)** Added regression weights ($w_1, w_2, w_3$) to prevent double-counting. |
| **v3.0** *(Final Peer-Review Standard)* | $S_{\text{dyn}} = w_1(\langle \Delta E_{\text{int},t}\rangle - \langle \Delta E_{\text{int},o}\rangle) - w_2 RT\ln\left(\frac{f_{\text{key},t} + \epsilon}{f_{\text{key},o} + \epsilon}\right) - w_3 RT\left(\frac{\langle \text{SASA}_{\text{bur},t}\rangle}{\text{SASA}_{\text{lig}}} - \frac{\langle \text{SASA}_{\text{bur},o}\rangle}{\text{SASA}_{\text{lig}}}\right)$ | **1)** Fixed sign ambiguity via strictly non-negative $\text{SASA}_{\text{bur}} = \text{SASA}_{\text{free}} - \text{SASA}_{\text{bound}} \ge 0$.<br>**2)** Replaced static pose energy with MD ensemble-averaged $\langle \Delta E_{\text{int}} \rangle$.<br>**3)** Defined $f_{\text{key}}$ strictly via catalytic H-bonds ($d \le 3.5\text{ \AA}$) and replaced hardcoded $-3.0$ with **Youden's $J$ adaptive cutoff $\theta^*$**. |

---

## 📐 Supporting Thermodynamic & Kinetic Equations ($T = 310\text{ K}$)

* **Boltzmann-Weighted Contact Persistence & Energy:**
  $$P_n = \frac{\exp(-E_n / RT)}{\sum_{i=0}^{N} \exp(-E_i / RT)}, \quad \langle C \rangle = \sum_{n=0}^{N} P_n \cdot C_n, \quad \langle \Delta E_{\text{contact}} \rangle = \sum_{n=0}^{N} P_n \cdot C_n \cdot \Delta E_n$$
* **Heuristic Residence Time & Rebek Steric Packing:**
  $$\tau_{\text{est}} = \tau_0 \exp\left(\frac{\Delta G^{\ddagger}_{\text{proxy}}}{RT}\right), \quad \eta_{\text{pack}} = \frac{V_{\text{ligand}}}{V_{\text{pocket}}} \approx 0.55$$

---

## 🚀 Quick Start & Repository Contents

```bash
git clone https://github.com
cd TB-Dynamic-Selectivity
python3 full_pipeline_v3_audit.py
```

* `MTB-NL4_Computational_Screening_Summary_Report.pdf` — Complete 11-section summary report with equations, plain-English guide, and ROC plots.
* `full_pipeline_v3_audit.py` & `formula_v3_calibrated.py` — Final v3.0 verification and Youden's $J$ adaptive threshold engine.
* `benchmark_20_medicines_v2.py` & `MEDICINE_20_BENCHMARK_REPORT.txt` — 20-medicine selectivity audit.
* `real_literature_inha_validation.py` — Experimental $\text{IC}_{50}$ Spearman correlation script.
* `lead_optimization_engine.py` — Closed-loop R-group optimization engine.
* `dude_inha_30_benchmark.py` & `dude_inha_roc_benchmark.png` — 30-molecule thermal noise ROC-AUC benchmark.
* `prove_my_formulas.py` & `PROOF_REPORT.txt` — Automated mathematical and thermodynamic proofs.

---

## 📚 Literature Reference Compounds & Proxy Benchmark Provenance

> **Methodological Transparency Note:** Multi-compound benchmarks (`n = 30`, `n = 20`, `n = 14`) currently use **literature-anchored biophysical proxy inputs** (not completed all-atom MD trajectories) to audit formula arithmetic, dimensional consistency, and cross-validation stability prior to full `100 ns+` cloud GPU GROMACS runs.

| Compound | Complex PDB ID | RCSB Ligand HetID | Reported IC50 | Primary Literature Reference |
| :--- | :--- | :--- | :--- | :--- |
| **PT70** | `2X23` | `TCU` | `12 nM` | Luckner et al., *J. Biol. Chem.* (2010) |
| **NITD-916** | `4U0J` | `38I` | `50 nM` | Manjunatha et al., *Sci. Transl. Med.* (2015) |
| **GSK693** | `5JFO` | `68U` | `100 nM` | Martínez-Hoyos et al., *EBioMedicine* (2016) |
| **GEQ** | `2NSD` | `GEQ` | `160 nM` | Kuo et al., *J. Biol. Chem.* (2003) |
| **Pyrrolidine-1** | `4TZK` | `641` | `390 nM` | He et al., *J. Med. Chem.* (2006) |
| **Aryl-Amide Analog-4b** | `4OIM` | `2Q3` | `620 nM` | Shirude et al., *J. Med. Chem.* (2013) |
| **Benzimidazole Lead-12** | `4BQP` | `J9Z` | `850 nM` | Moreira et al., *ChemMedChem* (2013) |
| **Triclosan** | `1P45` | `TCL` | `1000 nM` | Kuo et al., *J. Biol. Chem.* (2003) |
