import re
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

# -----------------------------
# чтение параметров из конфига
# -----------------------------

def read_grid_from_config(filename):

    Nx = None
    Ny = None

    with open(filename) as f:
        for line in f:

            if "N_x" in line:
                Nx = int(re.findall(r'\d+', line)[0])

            if "N_y" in line:
                Ny = int(re.findall(r'\d+', line)[0])

    return Nx, Ny


# -----------------------------
# чтение px py
# -----------------------------

def read_decomposition(filename):

    px = None
    py = None

    with open(filename) as f:
        for line in f:

            if "px" in line:
                px = int(re.findall(r'\d+', line)[0])

            if "py" in line:
                py = int(re.findall(r'\d+', line)[0])

    return px, py


# -----------------------------
# локальный размер
# -----------------------------

def local_size(N, coord, dim):

    base = N // dim
    rest = N % dim

    if coord < rest:
        return base + 1
    else:
        return base


# -----------------------------
# offset
# -----------------------------

def offset(N, coord, dim):

    base = N // dim
    rest = N % dim

    if coord < rest:
        return coord * (base + 1)
    else:
        return rest * (base + 1) + (coord - rest) * base


# -----------------------------
# визуализация
# -----------------------------

def plot_decomposition(Nx, Ny, px, py, save_path):

    fig, ax = plt.subplots(figsize=(6,6))

    rank = 0

    for ry in range(py):
        for rx in range(px):

            nx_loc = local_size(Nx, rx, px)
            ny_loc = local_size(Ny, ry, py)

            ox = offset(Nx, rx, px)
            oy = offset(Ny, ry, py)

            rect = patches.Rectangle(
                (ox, oy),
                nx_loc,
                ny_loc,
                linewidth=1,
                edgecolor="black",
                facecolor="lightblue",
                alpha=0.6
            )

            ax.add_patch(rect)

            ax.text(
                ox + nx_loc/2,
                oy + ny_loc/2,
                f"{rank}",
                ha="center",
                va="center"
            )

            rank += 1

    ax.set_xlim(0, Nx)
    ax.set_ylim(0, Ny)

    ax.set_xlabel("Nx")
    ax.set_ylabel("Ny")

    ax.set_title(f"Domain decomposition ({px} x {py})")

    ax.set_aspect("equal")

    plt.savefig(save_path, dpi=200)
    plt.close()
# -----------------------------
# main
# -----------------------------

if len(sys.argv) < 2:
    print("Usage: python3 picture.py <config_name>")
    sys.exit(1)

config_name = sys.argv[1]

config_path = Path("config_set") / f"{config_name}.toml"
decomp_path = Path("output") / config_name / "decomposition.txt"

Nx, Ny = read_grid_from_config(config_path)
px, py = read_decomposition(decomp_path)

save_path = Path("output") / config_name / "decomposition.png"

plot_decomposition(Nx, Ny, px, py, save_path)

print(f"Saved: {save_path}")