import os, numpy as np
from operator import itemgetter
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, FilterCatalog
from reportlab.lib.pagesizes import A4
from reportlab.platypus import KeepTogether, SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

smiles_str = "O=C(Nc1ccc(C(F)(F)F)cc1)c1cc(O)ccc1-c1ccccn1"
mol = Chem.MolFromSmiles(smiles_str)
formula, mw, clogp = rdMolDescriptors.CalcMolFormula(mol), Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
hbd, hba = rdMolDescriptors.CalcNumHBD(mol), rdMolDescriptors.CalcNumHBA(mol)
p = FilterCatalog.FilterCatalogParams()
p.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
pains_cnt = len(FilterCatalog.FilterCatalog(p).GetMatches(mol))
Draw.MolToFile(mol, "mtb_nl4_2d_structure.png", size=(420, 220))

RT, EPS = 0.001987204 * 310.0, 0.01
def S_v3(et, eo, ft, fo, st, so, sl, w1=1.0, w2=1.0, w3=1.5):
    return w1*(et - eo) - w2*RT*np.log((ft + EPS)/(fo + EPS)) - w3*RT*((st - so)/sl)

meds = [
    ("Isoniazid-NAD",1,-12.1,-7.2,0.88,0.05,440,240,520),("Ethionamide-NAD",1,-11.5,-7.0,0.81,0.07,420,255,515),
    ("Prothionamide-NAD",1,-11.4,-7.1,0.79,0.08,415,260,520),("PT70",1,-11.9,-7.4,0.85,0.05,430,250,510),
    ("NITD-916",1,-11.4,-7.1,0.78,0.06,415,260,520),("GSK693",1,-10.9,-6.8,0.71,0.07,405,265,515),
    ("Triclosan",1,-10.5,-7.2,0.68,0.12,395,290,480),("GEQ",1,-11.2,-7.0,0.75,0.08,410,270,510),
    ("Pyrrolidine-1",1,-10.8,-6.9,0.72,0.09,400,275,505),("AN12855",1,-11.1,-6.9,0.76,0.07,412,262,515),
    ("Paracetamol",0,-6.2,-6.4,0.12,0.18,210,230,350),("Ibuprofen",0,-7.5,-7.9,0.18,0.25,260,290,420),
    ("Atorvastatin",0,-9.4,-9.8,0.24,0.72,340,460,580),("Ketoconazole",0,-9.6,-9.5,0.27,0.68,350,450,560),
    ("Warfarin",0,-8.1,-8.3,0.20,0.32,280,310,450),("Metformin",0,-5.4,-5.5,0.08,0.10,180,190,310),
    ("Ciprofloxacin",0,-8.3,-8.1,0.21,0.24,290,305,470),("Fluconazole",0,-7.8,-7.9,0.17,0.22,270,285,440),
    ("Sorafenib",0,-9.8,-9.7,0.29,0.65,360,440,550),("Aspirin",0,-6.0,-6.1,0.10,0.14,195,210,340)
]
names = list(map(itemgetter(0), meds))
y_all = np.array(list(map(itemgetter(1), meds)))
s_all = np.array([S_v3(*m[2:]) for m in meds])

# Print Triclosan exact decomposition
tri_idx = names.index("Triclosan")
et, eo, ft, fo, st, so, sl = meds[tri_idx][2:]
t1, t2, t3 = (et - eo), -RT*np.log((ft+EPS)/(fo+EPS)), -1.5*RT*((st-so)/sl)
print(f"Triclosan Exact Score -> Term1: {t1:.3f}, Term2: {t2:.3f}, Term3: {t3:.3f} => Total S_dynamic = {s_all[tri_idx]:.3f} kcal/mol")
print(f"All 10 Positives sorted S_dynamic: {np.sort(s_all[y_all==1]).round(3)} (Least negative positive = {np.max(s_all[y_all==1]):.3f})")

pos_i, neg_i = np.where(y_all == 1).__getitem__(0), np.where(y_all == 0).__getitem__(0)
oof, thetas, missed = np.zeros_like(y_all), [], []
for k in range(5):
    te = np.concatenate([pos_i[k*2:(k+1)*2], neg_i[k*2:(k+1)*2]])
    tr = np.setdiff1d(np.arange(20), te)
    bj, bth = -1.0, 0.0
    for th in np.sort(np.unique(s_all[tr])):
        pr = (s_all[tr] <= th).astype(int)
        j = np.sum((pr==1)&(y_all[tr]==1))/np.sum(y_all[tr]==1) - np.sum((pr==1)&(y_all[tr]==0))/np.sum(y_all[tr]==0)
        if j > bj: bj, bth = j, th
    thetas.append(bth)
    oof[te] = (s_all[te] <= bth).astype(int)
    for i in te:
        if oof[i] != y_all[i]:
            missed.append(f"{names[i]} (S={s_all[i]:.3f} vs fold theta*={bth:.3f})")

n_correct, m_th = int(np.sum(oof == y_all)), np.mean(thetas)
miss_str = ", ".join(missed)
lit_14 = [
    (50.0,-11.4,-7.1,0.78,0.06,415,260,520),(12.0,-11.9,-7.4,0.85,0.05,430,250,510),(100.0,-10.9,-6.8,0.71,0.07,405,265,515),
    (1000.0,-10.5,-7.2,0.68,0.12,395,290,480),(390.0,-10.8,-6.9,0.72,0.09,400,275,505),(160.0,-11.2,-7.0,0.75,0.08,410,270,510),
    (620.0,-10.3,-6.8,0.66,0.11,390,280,495),(850.0,-10.1,-6.7,0.64,0.10,385,285,490),(45000.0,-8.2,-7.8,0.22,0.25,290,310,470),
    (38000.0,-9.6,-9.4,0.30,0.58,340,410,500),(60000.0,-7.9,-7.5,0.18,0.20,270,280,460),(85000.0,-7.4,-7.3,0.15,0.19,240,260,450),
    (32000.0,-9.2,-9.5,0.28,0.65,320,430,490),(52000.0,-8.5,-8.0,0.19,0.22,280,295,475)
]
pic50 = [-np.log10(x*1e-9) for x in map(itemgetter(0), lit_14)]
s14 = [S_v3(*x[1:]) for x in lit_14]
rho, pval = spearmanr(pic50, s14)
print(f"5-Fold CV: {n_correct}/20 ({n_correct*5.0:.1f}%) | Missed: {miss_str} | Spearman r={rho:.3f}, p={pval:.2e}")
pdf_name = "MTB-NL4_Computational_Screening_Summary_Report.pdf"
doc = SimpleDocTemplate(pdf_name, pagesize=A4, topMargin=26, bottomMargin=26)
styles = getSampleStyleSheet()
cell = ParagraphStyle("Cell", parent=styles["Normal"], fontSize=8.3, leading=10.8)
story = [
    Paragraph("<b>TB-Dynamic-Selectivity v3.0: Computational Screening &amp; Biophysical Audit Report</b>", styles["Title"]),
    Paragraph("<b>Author:</b> Navneet Dwivedi | <b>Target:</b> Mycobacterium tuberculosis InhA (PDB: 4TZK) | <b>Lead:</b> MTB-NL4", styles["Normal"]),
    Spacer(1, 8),
    Paragraph("<b>1. Executive Summary &amp; Benchmark Overview</b>", styles["Heading3"]),
    Paragraph(
        "Single-conformation docking (AutoDock Vina) ignores protein flexibility and the statistical weight of alternative poses. <b>TB-Dynamic-Selectivity</b> tests a physics-informed score (S_dynamic, v4.0) that combines an ensemble of docking poses, a smooth catalytic-contact term and a buried-surface term, relative to the InhA crystal ligand 641. <b>Current result (v4.0, docking poses only, no MD):</b> 25 InhA actives and 60 property- and similarity-matched decoys. On the structured benchmark (31 co-crystal actives vs 60 ZINC-linked decoys), static Vina AUC = 0.765 (95% CI 0.652&#8211;0.863) and the v4.0 score AUC = 0.767 (95% CI 0.658&#8211;0.865); permutation p &lt; 0.001 for both. Within MW-matched windows both methods reach AUC about 0.81&#8211;0.82. The two methods are not statistically distinguishable, so v4.0 is not shown to outperform static docking. Earlier in-sample proxy results (v1&#8211;v3 benchmarks, circular CV and Spearman checks) are no longer reported in this document and are not evidence of predictive performance.",
        styles["Normal"]
    ),
    Spacer(1, 6),
    Paragraph("<b>2. Formula Evolution Across 3 Iterations (v1.0 -&gt; v2.0 -&gt; v3.0 -&gt; v4.0)</b>", styles["Heading3"])
]
evo_data = [
    [Paragraph("<b>Version</b>", cell), Paragraph("<b>Core Structure</b>", cell), Paragraph("<b>Methodological Upgrades &amp; Fixes</b>", cell)],
    [Paragraph("<b>v1.0</b>", cell), Paragraph("ddG_dyn = (dG_t - dG_o) - RT*ln(f_t/f_o) - 1.5*RT*(dSASA_t/Vp_t - dSASA_o/Vp_o)", cell), Paragraph("Initial heuristic combining static Vina docking score, persistence, and burial (-3.0 cutoff).", cell)],
    [Paragraph("<b>v2.0</b>", cell), Paragraph("S_dyn = w1*(dG_t - dG_o) - w2*RT*ln((f_t+e)/(f_o+e)) - w3*RT*(dSASA_t - dSASA_o)/SASA_lig", cell), Paragraph("Fixed kcal/mol units (A^2/A^2), added epsilon=0.01 crash safety &amp; weights (w1, w2, w3).", cell)],
    [Paragraph("<b>v3.0</b>", cell), Paragraph("S_dyn = w1*(&lt;dE_t&gt; - &lt;dE_o&gt;) - w2*RT*ln((f_key,t+e)/(f_key,o+e)) - w3*RT*(&lt;SASA_bur,t&gt; - &lt;SASA_bur,o&gt;)/SASA_lig", cell), Paragraph("Resolved sign via SASA_bur &gt;= 0, MD ensemble &lt;dE_int&gt;, catalytic f_key (d &lt;= 3.5 A) &amp; Youden's J cutoff theta*.", cell)],
    [Paragraph("<b>v4.0</b>", cell), Paragraph("S = w1*(&lt;dE&gt;_t - &lt;dE&gt;_o) - w2*RT*ln((&lt;f_key&gt;_t+eps)/(&lt;f_key&gt;_o+eps)) - w3*RT*(&lt;SASA_bur&gt;_t/SASA_lig - &lt;SASA_bur&gt;_o/SASA_lig); tau_est = tau0*exp(-S/RT)", cell), Paragraph("Ensemble from 9 Vina poses (Boltzmann); smooth Tyr158 contact weight (1 at &lt;=3.0 A, 0 at &gt;=4.2 A) replaces hard 3.5 A step; eps = 0.02 chosen by principle; gate (dE_t &gt; 0 gives S = 1e6); reference = 641 (InhA ligand); tau_est heuristic, tracks S.", cell)]
]
t_evo = Table(evo_data, colWidths=(70, 200, 215))
t_evo.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#E2E8F0")),("GRID",(0,0),(-1,-1),0.6,colors.HexColor("#94A3B8")),("VALIGN",(0,0),(-1,-1),"TOP"),("PADDING",(0,0),(-1,-1),4)]))
story += [
    t_evo, Spacer(1, 6),
    Paragraph("<b>3. Core Mathematical Equations (v4.0) &amp; Conceptual Overview</b>", styles["Heading3"]),
    Image("formulas_v4_eqs.png", width=475, height=304),
    Spacer(1, 4)
]

story += [Paragraph("<b>3b. The three equations and what they mean</b>", styles["Heading3"]),
    Spacer(1, 2),
    Paragraph("<b>Equation 1: the score.</b> S compares the candidate with the reference molecule 641 (the crystal ligand of InhA). Three parts: (1) binding energy difference; (2) difference in how long the molecule stays near the catalytic residue Tyr158 (f_key, a value from 0 to 1); (3) difference in the fraction of its surface hidden in the protein (SASA_bur divided by the molecule free surface SASA_lig). S = 0 means equal to 641; negative means predicted to bind better than 641; positive means weaker.", styles["Normal"]), Spacer(1, 4),
    Paragraph("<i>Symbols:</i> t = candidate; ref = molecule 641. dE = Vina energy in kcal/mol (more negative is better). w1 = 1.0, w2 = 1.0, w3 = 1.5 are weights set by hand, not fitted. RT = 0.616 kcal/mol (thermal energy at 310 K). eps = 0.02 is a small number that keeps the logarithm defined when f_key is near zero; it was chosen from the smallest meaningful contact value, not from the benchmark. Gate: if the candidate average energy is above 0.0 kcal/mol, S is set to 1e6 (ranked last).", styles["Normal"]), Spacer(1, 6),
    Paragraph("<b>Equation 2: averaging over docking poses.</b> Docking gives 9 possible positions (poses) of one molecule, each with energy E_n. P_n turns the energies into probabilities, so the lowest-energy pose gets the largest weight and all weights add up to 1. Every quantity in Equation 1 (dE, f_key, SASA_bur) is first computed for each pose and then averaged with these weights.", styles["Normal"]), Spacer(1, 4),
    Paragraph("<i>Symbols:</i> E_n = Vina score of pose n. E_min = lowest score. X_n = the quantity for pose n. Sum over n runs over all 9 poses. Because RT is small, a pose 1.5 kcal/mol above the best gets about 9% of the best pose weight, so the best pose dominates the average.", styles["Normal"]), Spacer(1, 6),
    Paragraph("<b>Equation 3: residence-time index (heuristic).</b> tau_est grows as S becomes more negative. It is a rescaling of S, not a measured time: it has no time unit and has not been checked against experiment. Because it is a monotonic function of S, it ranks molecules in exactly the same order as S, so it adds no new information and its AUC equals that of S.", styles["Normal"]), Spacer(1, 6),
    Paragraph("<b>Benchmark (same system).</b> 25 InhA actives (6 verified, 19 from ChEMBL, IC50 at most 1 uM; CHEMBL3125101 excluded because Vina does not support boron) and 60 DUD-E decoys (property-matched; Morgan Tanimoto at most 0.23 to every active). On the structured benchmark (31 actives, 60 ZINC-linked decoys) Vina AUC = 0.765 (95% CI 0.652-0.863) and S_v4 AUC = 0.767 (95% CI 0.658-0.865); permutation p below 0.001 for both. The earlier eps sensitivity (0.001, 0.02, 0.05) was computed on the superseded 25-active set and has not been recomputed. The two methods are not statistically different, so v4.0 is not shown to beat plain docking.", styles["Normal"]), Spacer(1, 6),
    Paragraph("<b>Limits.</b> All averages come from 9 docking poses, not molecular dynamics, so the 310 K statement applies only to docking poses. The reference is one molecule (641), not a human off-target, so S measures similarity to 641, not selectivity. Weights and eps are set by hand. The burial term is about 5% of S. The smooth contact weight replaced an earlier hard 3.5 A cutoff that caused a jump of about 4 kcal/mol for a 0.07 A change. At RT = 0.616, one pose gets most of the weight, so the ensemble is effectively one pose.", styles["Normal"]), Spacer(1, 4)]

secs = [
          ("4. Simulation Protocol &amp; Benchmark Data Provenance (Transparency Disclosure)", "<b>Target System &amp; MD Setup:</b> M. tuberculosis InhA crystal structure (PDB: 4TZK, chain A, co-crystallized with NAD+ and pyrrolidine carboxamide inhibitor 641; He et al., J. Med. Chem. 2006), parameterized with AMBER99SB-ILDN (protein) and GAFF2/AM1-BCC (ligand). Short GROMACS runs (minimization, 100 ps NVT at 310 K, 400 ps NPT at 310 K) were performed for three InhA-NAD-ligand complexes; only final structures were kept, no trajectory analysis has been performed, and no result in this report is derived from these runs. Production sampling (100 ns or longer) remains planned.<br/><b>Benchmark (v4.0):</b> 25 InhA actives (6 structurally verified reference compounds plus 19 ChEMBL-curated M. tuberculosis inhibitors with IC50 at most 1 uM; CHEMBL3125101 excluded because Vina does not support boron) and 60 DUD-E decoys (property-matched; Morgan Tanimoto at most 0.23 to every active). Potency values are summarized in Table 2 with their audit status. Activity values were used only as ordering proxies and have not all been re-verified against full-text tables."),
("__TABLE2__", ""),
    ("Docking validation: redocking control (pose reproduction)", "Redocking of verified crystal ligands into their own receptors (Vina, exhaustiveness 32, 9 modes, box = crystal ligand extent + 8 A) was used as a pose-reproduction control. Heavy-atom RMSD was computed with a graph-based, symmetry-aware atom mapping against the crystal ligand, without alignment. PT119 (4OIM): 0.65 A; PT70 (2X23): 0.71 A; 641 (4TZK, common receptor): 0.23 A (reproduced). Triclosan (1P45, residue 400, chain A): 5.22 A; NITD-916 (4R9S): 5.92 A (centroid 2.5 A, orientation not reproduced); GSK625 (5JFO): 8.89 A; 4PI (2NSD): 9.27 A (pose not reproduced). Conclusion: docking reproduces the crystal pose for three of seven verified ligands and places four others in the pocket with altered orientation or position. Docking poses are therefore treated as ranking proxies, not reliable binding modes, for the full set. Triclosan was docked from residue 400 only; chain A also contains residue 450 and chain B a third copy, which were excluded."),
("5. Top Screened Lead Candidate: MTB-NL4 (RDKit Audit &amp; 2D Structure)",
     f"<b>Candidate ID:</b> MTB-NL4 (Novel_Lead_4_5CF3) | <b>Target:</b> M. tuberculosis InhA (PDB: 4TZK)<br/>"
     f"<b>IUPAC Name:</b> 5-hydroxy-2-(pyridin-2-yl)-N-[4-(trifluoromethyl)phenyl]benzamide<br/>"
     f"<b>Canonical SMILES:</b> {smiles_str}<br/>"
     f"<b>RDKit Physicochemical Audit:</b> Formula: {formula} | MW: {mw:.2f} g/mol | cLogP: {clogp:.2f} (RDKit Crippen) | HBD: {hbd}, HBA: {hba} (RDKit rdMolDescriptors; Lipinski compliant)<br/>"
     f"<b>Structural Alerts:</b> {pains_cnt} PAINS alerts (RDKit FilterCatalog PAINS) | <b>AutoDock Vina Mode 1 Score:</b> −8.60 kcal/mol (receptor without NAD, inha_receptor.pdbqt; the same ligand gave −8.56 kcal/mol on the NAD-containing receptor inha_clean.pdb). The earlier −11.703 kcal/mol run used an unverified setup and is not used."),
    ("6. Benchmark summary (v4.0, structured actives vs ZINC-linked decoys)", "Benchmark: 31 structured InhA co-crystal inhibitors (wild-type InhA, RCSB-verified; cofactor and ions removed; EEH excluded for stereo-parsing reasons) and 60 DUD-E decoys with ZINC IDs linked to their own Vina and v4.0 scores (decoy_scores_by_zinc.tsv). Full set: Vina AUC 0.765 (95% CI 0.652-0.863); v4.0 S AUC 0.767 (95% CI 0.658-0.865); permutation p < 0.001 for both. MW-matched subsets: MW 280-420 (26 actives, 57 decoys) Vina 0.812, S 0.823; MW 260-440 (28 actives, 60 decoys) Vina 0.819, S 0.817; p < 0.001 in all windows. The two methods are indistinguishable in every window (overlapping CIs), so v4.0 is not shown to outperform static docking. Earlier results (25 actives, Vina 0.858, S 0.822) were obtained with a decoy index that was not verified against ZINC IDs and are superseded. Sensitivity to eps was reported on the earlier set and has not been recomputed on the structured benchmark."),
("__ROC__", ""),
    ("7. Analogue design (computational property screen only)", "Five analogues of MTB-NL4 were screened with RDKit for physicochemical properties and PAINS/Brenk alerts. These are design proposals only: they were not docked, synthesized, or tested, and no activity or selectivity is claimed for them."),
    ("8. Open-Source Repository &amp; Reproducibility", "All source code, v1-v3 evolution scripts, and benchmark reports are publicly hosted at: https://github.com/navneetdwivedi147-code/TB-Dynamic-Selectivity"),
    ("9. Planned work (not yet performed)", "Production molecular dynamics (100 ns or longer, all-atom GROMACS, NPT, 310 K) for the InhA compounds in Table 2 and the benchmark set is planned, together with nested cross-validation for the score weights and a preprint. No result from these planned steps is reported here. Computing resources are a constraint on this plan; dedicated GPU time has not been secured."),
    ("Appendix A. beta-1 AR (7BVQ) pilot (separate project; not part of the InhA result)", "<b>Setup:</b> 7BVQ (beta-1 AR, carazolol-bound) with ND-2606-OPT complex (SMILES: CC(C)[NH2+]C[C@H](O)COc1cccc2c1c(=O)n(CCO)cc2; identity, origin and novelty unknown to the authors; no prior-art claim is made), TIP3P solvent only, 0.25 ns per replicate (251 frames).<br/><b>Ligand RMSD (backbone fit):</b> Replicate 1 mean 1.75 A (max 2.66 A); Replicate 2 mean 1.33 A (max 2.34 A).<br/><b>Asp1138 contact (ligand N/O3 to Asp1138 carboxylate):</b> mean minimum distance 4.14 A (Rep 1) and 4.28 A (Rep 2); frames with distance &lt; 3.5 A: 0.0% in both runs; frames with distance &lt;= 4.0 A: 27.1% (Rep 1) and 6.4% (Rep 2). The docking contact (2.82 to 2.85 A) was not maintained in either run.<br/><b>Protein C-alpha RMSF:</b> mean 0.97 A (Rep 1) and 1.06 A (Rep 2); Asp1138 0.64 A (Rep 1) and 0.49 A (Rep 2); pocket residues 1136-1140 show low flexibility.<br/><b>Interpretation:</b> The ionic contact did not form in either replicate, but the close-contact fraction differs between runs, so this result is not reproducible at 0.25 ns."),
("__REFS__", ""),
]
for t, b in secs:
    if t == "__ROC__":
        import os
        if os.path.exists("benchmark_roc_v4.png"):
            story += [Spacer(1, 4), Image("benchmark_roc_v4.png", width=300, height=266), Paragraph("<b>Figure 2: ROC curves for static Vina and the v4.0 score (31 co-crystal actives vs 60 ZINC-linked decoys; Vina AUC 0.765, v4.0 AUC 0.767).</b>", styles["Normal"]), Spacer(1, 6)]
        continue
    if t == "__TABLE2__":
        from sec4_table import build_table2
        story += build_table2(styles)
        continue
    if t == "__REFS__":
        from sec4_table import build_refs
        story += build_refs(styles)
        continue
    story += [Paragraph(f"<b>{t}</b>", styles["Heading3"]), Paragraph(b, styles["Normal"]), Spacer(1, 4)]
if os.path.exists("mtb_nl4_2d_structure.png"):
    story += [Spacer(1, 4), Paragraph("<b>Figure 1: 2D Chemical Structure of Lead Candidate MTB-NL4 (Novel_Lead_4_5CF3)</b>", styles["Normal"]), Image("mtb_nl4_2d_structure.png", width=250, height=135)]
doc.build(story)
print("Generated MTB-NL4_Computational_Screening_Summary_Report.pdf")

