import numpy as np

RT = 0.001987204 * 310.0  # 0.616 kcal/mol at 310 K

# 5 MD States: Total State Energy (E), Key-Residue Contact (C: 1 or 0), Local Contact Energy (dE_local)
E        = np.array([-11.8, -11.5, -10.2, -9.5, -8.8])
C        = np.array([1, 1, 0, 1, 0])
dE_local = np.array([-4.2,  -3.9,   0.0, -2.1,  0.0])

# Boltzmann Probabilities
weights = np.exp(-(E - E.min()) / RT)
P = weights / np.sum(weights)

# Metric 1: Boltzmann-Weighted Contact Persistence (<C> %)
boltzmann_contact_pct = np.sum(P * C) * 100.0

# Metric 2: Boltzmann-Weighted Contact Energy (<dE_contact> kcal/mol)
boltzmann_contact_energy = np.sum(P * C * dE_local)

print("\n=== 🔬 DUAL BOLTZMANN CONTACT & ENERGY AUDIT ===")
print(f"1. Boltzmann-Weighted Contact Persistence (<C>)   : {boltzmann_contact_pct:.1f} %")
print(f"2. Boltzmann-Weighted Contact Energy (<dE_contact>): {boltzmann_contact_energy:.3f} kcal/mol")
print("===================================================\n")
