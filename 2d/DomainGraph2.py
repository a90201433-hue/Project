import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
import numpy as np
import sys

config = sys.argv[1]

file = Path("output") / config / "domains.txt"

domains = []

with open(file) as f:
    next(f)
    for line in f:
        vals = list(map(float, line.split()))
        domains.append(vals)

p = len(domains)
colors = plt.cm.tab20(np.linspace(0,1,p))

fig, ax = plt.subplots(figsize=(7,7))

for i, d in enumerate(domains):

    rank = int(d[0])

    xmin, xmax, ymin, ymax = d[3], d[4], d[5], d[6]
    xmin_g, xmax_g, ymin_g, ymax_g = d[7], d[8], d[9], d[10]

    color = colors[i]

    # физическая область
    rect = patches.Rectangle(
        (xmin, ymin),
        xmax - xmin,
        ymax - ymin,
        facecolor=color,
        edgecolor="black",
        alpha=0.6,
        linewidth=1.5
    )

    # ghost область
    rect_g = patches.Rectangle(
        (xmin_g, ymin_g),
        xmax_g - xmin_g,
        ymax_g - ymin_g,
        fill=False,
        edgecolor=color,
        linestyle="--",
        linewidth=2
    )

    ax.add_patch(rect_g)
    ax.add_patch(rect)

    ax.text(
        (xmin + xmax)/2,
        (ymin + ymax)/2,
        str(rank),
        ha="center",
        va="center",
        fontsize=10
    )

ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.title("Domain decomposition with ghost cells")

plt.savefig(Path("output") / config / "domains_ghost.png", dpi=200)