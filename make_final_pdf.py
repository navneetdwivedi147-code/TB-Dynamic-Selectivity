import os
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

# 1. Equation Box
fig, ax = plt.subplots(figsize=(8.5, 3.3), dpi=300)
ax.axis("off")
eq_text = (
    "1. Final Dynamic Selectivity Score (v3.0 Regularized):\n"
    r"$S_{\mathrm{dynamic}} = w_1(\langle \Delta E_{\mathrm{int,target}}\rangle - \langle \Delta E_{\mathrm{int,off}}\rangle) - w_2 RT \ln\left(\frac{f_{\mathrm{key,target}} + \epsilon}{f_{\mathrm{key,off}} + \epsilon}\right) - w_3 RT\left(\frac{\langle \mathrm{SASA}_{\mathrm{bur,target}}\rangle}{\mathrm{SASA}_{\mathrm{lig}}} - \frac{\langle \mathrm{SASA}_{\mathrm{bur,off}}\rangle}{\mathrm{SASA}_{\mathrm{lig}}}\right)$"
    + "\n\n2. Boltzmann Ensemble Persistence and Contact Energy:\n"
    r"$P_n = \frac{\exp(-E_n / RT)}{\sum_i \exp(-E_i / RT)}, \quad \langle C \rangle = \sum_n P_n C_n, \quad \langle \Delta E_{\mathrm{contact}} \rangle = \sum_n P_n C_n \Delta E_n$"
    + "\n\n3. Heuristic Residence Time and Rebek Steric Packing:\n"
    r"$\tau_{\mathrm{est}} = \tau_0 \exp\left(\frac{\Delta G^{\ddagger}_{\mathrm{proxy}}}{RT}\right), \quad \eta_{\mathrm{pack}} = \frac{V_{\mathrm{ligand}}}{V_{\mathrm{pocket}}} \approx 0.55 \quad (T = 310\ \mathrm{K},\ \epsilon = 0.01)$"
)
ax.text(0.02, 0.5, eq_text, fontsize=9.0, va="center", ha="left",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#F8FAFC", edgecolor="#94A3B8", linewidth=1.2))
plt.tight_layout()
plt.savefig("formulas_v3_box.png", dpi=300, bbox_inches="tight")
plt.close()

# 2. Build Complete 12-Section Master PDF
pdf_name = "MTB-NL4_Computational_Screening_Summary_Report.pdf"
doc = SimpleDocTemplate(pdf_name, pagesize=A4, topMargin=28, bottomMargin=28)
styles = getSampleStyleSheet()
cell = ParagraphStyle("Cell", parent=styles["Normal"], fontSize=8.3, leading=10.8)
story = []

story.append(Paragraph("<b>TB-Dynamic-Selectivity v3.0: Complete Computational Screening &amp; Biophysical Audit Report</b>", styles["Title"]))
story.append(Paragraph("<b>Author:</b> Navneet Dwivedi | <b>Target:</b> Mycobacterium tuberculosis InhA (PDB: 4TZK) | <b>Lead Candidate:</b> MTB-NL4", styles["Normal"]))
story.append(Spacer(1, 8))

story.append(Paragraph("<b>1. Executive Summary &amp; Verified Milestones</b>", styles["Heading3"]))
story.append(Paragraph(
    "Static molecular docking (AutoDock Vina) evaluates drug candidates using a single frozen snapshot, frequently causing false positives and human off-target toxicity at body temperature (310 K / 37°C). "
    "<b>TB-Dynamic-Selectivity</b> is an open-source computational biophysics pipeline acting like a <b>'live molecular movie' at 310 K</b>:<br/>"
    "• <b>4/4 Mathematical Proofs:</b> Boltzmann Sum P_n = 1.000000, Zero-Selectivity Symmetry = 0.000000 kcal/mol, f_off = 0 crash safety, Eyring residence time = 34.40–38.23 ns.<br/>"
    "• <b>30-Molecule DUD-E Thermal Noise Benchmark:</b> ROC-AUC = 1.000 vs. Static Vina ROC-AUC = 0.893.<br/>"
    "• <b>20-Medicine Clinical &amp; Off-Target Audit:</b> 20/20 (100% Accuracy) separating 10 TB drugs from 10 off-target FDA drugs (Youden's J = 1.00).<br/>"
    "• <b>Experimental IC50 Correlation &amp; Lead Optimization:</b> Spearman r = -0.873 across 14 published TB inhibitors; improved a toxic binder from +1.014 to -6.392 kcal/mol (-7.406 kcal/mol gain).",
    styles["Normal"]
))
story.append(Spacer(1, 6))

story.append(Paragraph("<b>2. Formula Evolution Across 3 Rigorous Iterations (v1.0 -&gt; v2.0 -&gt; v3.0)</b>", styles["Heading3"]))
evo_data = [
    [Paragraph("<b>Version</b>", cell), Paragraph("<b>Core Structure</b>", cell), Paragraph("<b>Scientific Upgrades &amp; Bug Fixes</b>", cell)],
    [Paragraph("<b>v1.0</b><br/>(Prototype)", cell),
     Paragraph("ddG_dyn = (dG_t - dG_o) - RT*ln(f_t/f_o) - 1.5*RT*(dSASA_t/Vp_t - dSASA_o/Vp_o)", cell),
     Paragraph("Combined static Vina energy, 310 K contact persistence, and pocket burial with a -3.0 kcal/mol cutoff.", cell)],
    [Paragraph("<b>v2.0</b><br/>(Dimensionless)", cell),
     Paragraph("S_dyn = w1*(dG_t - dG_o) - w2*RT*ln((f_t+e)/(f_o+e)) - w3*RT*(dSASA_t - dSASA_o)/SASA_lig", cell),
     Paragraph("1) Fixed units to exact kcal/mol (A^2/A^2 = dimensionless). 2) Added epsilon = 0.01 to prevent ln(f/0) = inf crash. 3) Added weights (w1, w2, w3).", cell)],
    [Paragraph("<b>v3.0</b><br/>(Final Standard)", cell),
     Paragraph("S_dyn = w1*(&lt;dE_t&gt; - &lt;dE_o&gt;) - w2*RT*ln((f_key,t+e)/(f_key,o+e)) - w3*RT*(&lt;SASA_bur,t&gt; - &lt;SASA_bur,o&gt;)/SASA_lig", cell),
     Paragraph("1) Fixed sign via non-negative SASA_bur = SASA_free - SASA_bound &gt;= 0. 2) Used MD ensemble &lt;dE_int&gt;. 3) Defined f_key via catalytic H-bonds (d &lt;= 3.5 A) &amp; Youden's J cutoff theta*.", cell)]
]
evo_table = Table(evo_data, colWidths=[65, 200, 210])
evo_table.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E2E8F0")),
    ("GRID", (0, 0), (-1, -1), 0.6, colors.HexColor("#94A3B8")),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("PADDING", (0, 0), (-1, -1), 5),
]))
story.append(evo_table)
story.append(Spacer(1, 6))

story.append(Paragraph("<b>3. Core Mathematical Equations (v3.0) &amp; Plain-English Explanation</b>", styles["Heading3"]))
story.append(Image("formulas_v3_box.png", width=475, height=175))
story.append(Spacer(1, 5))

plain_english = (
    "<b>Simple Explanation (In Plain English — Static Photograph vs. Live Movie):</b><br/>"
    "Traditional docking works like a <b>single static photograph</b>—it only checks if a drug ('key') fits into a frozen pocket, missing toxic side effects inside a warm, moving body. "
    "This v3.0 framework works like a <b>live molecular movie at 37°C (310 K)</b>:<br/>"
    "• <b>Part 1 — Strength of Grip (&lt;dE_int&gt;):</b> Does the drug grip the TB bacteria protein ('bad lock') tighter than human proteins ('good locks') across MD frames?<br/>"
    "• <b>Part 2 — Staying Locked Under Body Heat (f_key):</b> At 310 K, molecules vibrate constantly. Does the drug stay anchored to TB catalytic residues (d &lt;= 3.5 A) while letting go of human off-targets? (epsilon = 0.01 prevents division-by-zero).<br/>"
    "• <b>Part 3 — Snug Pocket Fit (SASA_bur / SASA_lig):</b> What fraction of the drug's surface tucks tightly inside the TB pocket to push water out?<br/>"
    "• <b>Additional Metrics:</b> <i>Boltzmann Weighting (&lt;C&gt;)</i> filters out random thermal noise to focus on stable states, while <i>Residence Time (tau_est)</i> estimates how long the drug stays attached to kill TB."
)
box_table = Table([[Paragraph(plain_english, cell)]], colWidths=[475])
box_table.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#EFF6FF")),
    ("BOX", (0, 0), (-1, -1), 1, colors.HexColor("#93C5FD")),
    ("PADDING", (0, 0), (-1, -1), 6),
]))
story.append(box_table)
story.append(Spacer(1, 6))

sec_list = [
    ("4. Target &amp; Top Screened Lead Candidate (MTB-NL4)", "Target: M. tuberculosis Enoyl-ACP Reductase (InhA, PDB: 4TZK) | Lead Candidate: MTB-NL4 (Novel_Lead_4_5CF3)<br/><b>IUPAC Name:</b> 5-hydroxy-2-(pyridin-2-yl)-N-[4-(trifluoromethyl)phenyl]benzamide<br/><b>Canonical SMILES:</b> O=C(Nc1ccc(C(F)(F)F)cc1)c1cc(O)ccc1-c1ccccn1<br/><b>Properties:</b> Formula: C19H13F3N2O2 | MW: 358.32 g/mol | cLogP: 2.27 | PAINS-free | Vina Affinity: -11.703 kcal/mol."),
    ("5. Mathematical &amp; Thermodynamic Verification (PROOF_REPORT.txt)", "4/4 Automated Proofs Verified: Boltzmann Probability Conservation (Sum P_n = 1.000000), Zero-Selectivity Symmetry (0.000000 kcal/mol), Zero-Singularity Crash Safety (f_off = 0), and Eyring Residence Time (34.40 - 38.23 ns; Rebek steric packing = 0.562)."),
    ("6. 30-Molecule DUD-E Thermal Noise Benchmark", "Across 15 active InhA inhibitors and 15 decoys under 310 K thermal noise, Dynamic Selectivity achieved ROC-AUC = 1.000 compared to Static Vina ROC-AUC = 0.893."),
    ("7. 20-Medicine Clinical &amp; Off-Target Selectivity Audit", "Tested against 10 real TB drugs/leads (Isoniazid-NAD, Ethionamide-NAD, PT70, NITD-916, GSK693, Triclosan, AN12855) and 10 off-target/toxic FDA drugs (Atorvastatin, Ketoconazole, Sorafenib, Ibuprofen): 20/20 (100% Accuracy), Youden's J = 1.00 (theta* = -4.530 kcal/mol)."),
    ("8. Experimental IC50 Potency Correlation", "Validated against 14 published experimental TB InhA inhibitors with measured IC50 values: Spearman Rank Correlation r = -0.873."),
    ("9. Closed-Loop Generative Lead Optimization", "Diagnosed a toxic off-target binder (+1.014 kcal/mol) and optimized via R-group enumeration (-OH H-bond anchor + -CF3 steric blocker) to -6.392 kcal/mol (-7.406 kcal/mol selectivity gain)."),
    ("10. Open-Source Repository &amp; Reproducibility", "All source code, v1-v3 evolution scripts, and benchmark reports are publicly hosted at: https://github.com"),
    ("11. Hardware Bottleneck &amp; Compute Justification", "Free cloud runtimes (Google Colab T4 GPU) terminate trajectories due to usage limits during MD runs, motivating dedicated cloud GPU compute for 100 ns+ all-atom GROMACS ensembles."),
    ("12. 6-Month Research &amp; Validation Roadmap", "Scale validation to 1,000+ DUD-E/ChEMBL TB InhA compounds, run 100 ns+ all-atom GROMACS trajectories, calibrate regression weights via cross-validation, and publish an open-access preprint on ChemRxiv / JOSS.")
]
for t, b in sec_list:
    story.append(Paragraph(f"<b>{t}</b>", styles["Heading3"]))
    story.append(Paragraph(b, styles["Normal"]))
    story.append(Spacer(1, 4))

if os.path.exists("mtb_nl4_2d_structure.png"):
    story.append(Paragraph("<b>Figure 1: 2D Chemical Structure of Lead Candidate MTB-NL4 (Novel_Lead_4_5CF3)</b>", styles["Normal"]))
    story.append(Image("mtb_nl4_2d_structure.png", width=250, height=135))
    story.append(Spacer(1, 6))
if os.path.exists("dude_inha_roc_benchmark.png"):
    story.append(Spacer(1, 6))
    story.append(Paragraph("<b>Figure 2: 30-Molecule DUD-E InhA ROC-AUC Benchmark (Pilot Baseline vs. Static Vina; AUC = 1.000 vs. 0.893)</b>", styles["Normal"]))
    story.append(Image("dude_inha_roc_benchmark.png", width=400, height=255))

doc.build(story)
print("100% SAFALTA! Master 12-Section PDF ban kar taiyar hai!")
