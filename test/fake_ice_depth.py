#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

depth = [5, 7.5, 10, 15, 20, 30, 40, 50, 60]
events = [377, 357, 336, 309, 279, 205, 144, 100, 80]

# create our plot
fig, ax = plt.subplots(figsize=(4, 4))

colors = plt.cm.get_cmap("Blues_r")(np.linspace(0, 1, 3))

ax.plot(depth, events, color=colors[0, :], ls="solid")
ax.scatter(depth, events, color=colors[0, :], s=5.2)

# add some labels
ax.set(
    xlabel="Ice Depth [m]",
    ylabel="Number of Events / 2yr",
    xticks=[0, 10, 20, 30, 40, 50, 60],
    ylim=[0, 400],
)

ax.axhline(275, ls="dashed", c="darkgrey", label="Requirement", zorder=-3)

# force the bottom to be zero
# ax.set_ylim(bottom=0.)

# and set a title
plt.title("Ice Depth")

# add our grid and legend
ax.set_axisbelow(True)
ax.grid(which="both", linestyle="dashed")

ax.legend(loc="upper right")

fig.savefig("ice_depth_scaled.png")
