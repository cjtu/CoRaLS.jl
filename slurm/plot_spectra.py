#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import ast

# ─── Read in all the outputs ────────────────────────────────────────────────
reflected_counts = []
reflected_errs   = []
direct_counts    = []
direct_errs      = []
altitudes        = []

for i in range(1, 21):
    altitudes.append(5 * i)                                # 5, 10, …, 100 km
    fname = f"out/old/corals_{i}km_10m.out"
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("r_spectra"):
                reflected_counts.append(
                    ast.literal_eval(line[12:].strip())
                )
            elif line.startswith("r_error"):
                reflected_errs.append(
                    ast.literal_eval(line[12:].strip())
                )
            elif line.startswith("d_spectra"):
                direct_counts.append(
                    ast.literal_eval(line[12:].strip())
                )
            elif line.startswith("d_error"):
                direct_errs.append(
                    ast.literal_eval(line[12:].strip())
                )

# ─── Convert to numpy and transpose so shape = (nE, nAlt) ───────────────────
rc_arr = np.array(reflected_counts)  # shape (nAlt, nE)
re_arr = np.array(reflected_errs)
dc_arr = np.array(direct_counts)
de_arr = np.array(direct_errs)

rc_transpose = rc_arr.T  # now (nE, nAlt)
re_transpose = re_arr.T
dc_transpose = dc_arr.T
de_transpose = de_arr.T

# ─── Define your energy bins (5 bins, 0.5 → 500 EeV linear) ─────────────────
nE = rc_transpose.shape[0]  # should be 5
Energies = [0.5 + 124.875 * j for j in range(nE)]  # [0.5, 125.375, 250.25, 375.125, 500.0]

# ─── Create a yellow→purple gradient ────────────────────────────────────────
cmap = mcolors.LinearSegmentedColormap.from_list("my_gradient", ["orange", "purple"], N=256)
colors = [cmap(j/(nE-1)) for j in range(nE)]

# ─── Plot reflected spectra vs altitude ─────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
for j in range(nE):
    ax.errorbar(
        altitudes,                  # x = altitudes
        rc_transpose[j],            # y = reflected counts at energy j
        yerr=re_transpose[j],       # errorbars
        color=colors[j],
        label=f"{Energies[j]:.3f} EeV",
        marker='o',
        linestyle='-'
    )
ax.set_xlabel("Altitude (km)")
ax.set_ylabel("Reflected Rate (# / yr)")
ax.set_title("Reflected CR Rate vs Altitude (10 m ice depth)")
ax.legend(title="Energy (EeV)")
plt.tight_layout()
plt.savefig("reflected_counts_10km.png", bbox_inches='tight')
# ─── Plot direct spectra vs altitude ────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 6))
for j in range(nE):
    ax2.errorbar(
        altitudes,
        dc_transpose[j],
        yerr=de_transpose[j],
        color=colors[j],
        label=f"{Energies[j]:.3f} EeV",
        marker='o',
        linestyle='-'
    )
ax2.set_xlabel("Altitude (km)")
ax2.set_ylabel("Direct Rate")
ax2.set_title("Direct Spectra vs Altitude")
ax2.legend(title="Energy (EeV)")
plt.tight_layout()

#plt.show()

