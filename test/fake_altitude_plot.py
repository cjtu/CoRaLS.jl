#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

altitudes = [20, 30, 40, 60, 80]
N8 = [308, 288, 282, 272, 266]
N12 = [437, 421, 409, 398, 387]
N16 = [639, 600, 578, 550, 531]
events = [N8, N12, N16]
Nants = [8, 12, 16]

# create our plot
fig, ax = plt.subplots(figsize=(4, 4))

colors = plt.cm.get_cmap("Blues_r")(np.linspace(0, 1, len(Nants) + 1))

for i in [0, 1, 2]:
    ax.plot(altitudes, events[i], color=colors[i, :], ls="solid", label=f"N={Nants[i]}")
    ax.scatter(altitudes, events[i], color=colors[i, :], s=5.2)

# add some labels
ax.set(
    xlabel="Orbit Altitude [km]",
    ylabel="Number of Events / 2yr",
    xticks=altitudes,
    ylim=[250, 700],
)

ax.axhline(275, ls="dashed", c="darkgrey", label="Requirement", zorder=-3)

# force the bottom to be zero
# ax.set_ylim(bottom=0.)

# and set a title
# plt.title("1 GHz BW @ 20km w/ Ice Roughness")

# add our grid and legend
ax.set_axisbelow(True)
ax.grid(which="both", linestyle="dashed")

ax.legend(loc="upper right")

fig.savefig("altitude_comparison_scaled.png")
