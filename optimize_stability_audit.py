from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog, QED

# Initialize both PAINS and Brenk catalogs
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog.FilterCatalog(params)

designs = [
    ("MTB-NL4 (Current Lead)",
     "O=C(Nc1ccc(C(F)(F)F)cc1)c1cc(O)ccc1-c1ccccn1",
     "Baseline: High cLogP (4.73) & exposed phenolic -OH"),

    ("Opt-1: Replace -CF3 with -CN",
     "O=C(Nc1ccc(C#N)cc1)c1cc(O)ccc1-c1ccccn1",
     "Removes -CF3, adds -CN -> cuts lipophilicity & hERG risk"),

    ("Opt-2: Pyridine -> Pyrimidine",
     "O=C(Nc1ccc(C(F)(F)F)cc1)c1cc(O)ccc1-c1cncnc1",
     "Adds 1 ring N -> lowers cLogP while keeping -CF3"),

    ("Opt-3: Dual (-CN + Pyrimidine + F)",
     "O=C(Nc1ccc(C#N)cc1F)c1cc(O)ccc1-c1cncnc1",
     "Property-filtered proposal: cLogP ~3.1, ortho-F blocks oxidation"),

    ("Opt-4: Replace -OH with -NHSO2Me",
     "CS(=O)(=O)Nc1ccc(-c2cncnc2)c(C(=O)Nc2ccc(C#N)cc2)c1",
     "Removes phenolic -OH to prevent rapid liver glucuronidation"),
]

print(f"{'Candidate':<33} | {'MW':>5} | {'cLogP':>5} | {'RotB':>4} | {'PAINS/Brenk':>11} | {'QED':>4} | {'Stability Verdict'}")
print("-" * 105)

for name, smi, note in designs:
    mol = Chem.MolFromSmiles(smi)
    mw = Descriptors.MolWt(mol)
    clogp = Descriptors.MolLogP(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    alerts = len(catalog.GetMatches(mol))
    qed = QED.qed(mol)

    if 2.0 <= clogp <= 3.6 and alerts == 0 and mw <= 400 and rotb <= 6:
        verdict = "passes simple property filters"
    elif clogp > 4.2:
        verdict = "HIGH cLogP RISK"
    else:
        verdict = "MODERATE"

    print(f"{name:<33} | {mw:5.1f} | {clogp:5.2f} | {rotb:4d} | {alerts:11d} | {qed:4.2f} | {verdict}")
