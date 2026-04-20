import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import glob

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "figure.dpi": 120,
})

# Геометрия: область 7×4 см, канал 2×1 см снизу слева
X_MIN, X_MAX = 0.0, 7.0
Y_MIN, Y_MAX = 0.0, 4.0
X_CORNER = 2.0
Y_CORNER = 1.0
WALL_RHO = 100.0
P_DET    = 0.3562
D_CJ     = 0.88
X_DET    = 0.15
P_THRESH = P_DET * 0.05

ASPECT = (X_MAX - X_MIN) / (Y_MAX - Y_MIN)


def read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"Ошибка чтения {path}: {e}")
        return None


def to_grid(df, field):
    xs = np.sort(df["x"].unique())
    ys = np.sort(df["y"].unique())
    grid = np.full((len(ys), len(xs)), np.nan)
    xi_idx = {v: i for i, v in enumerate(xs)}
    yi_idx = {v: i for i, v in enumerate(ys)}
    for _, row in df.iterrows():
        grid[yi_idx[row["y"]], xi_idx[row["x"]]] = row[field]
    return xs, ys, grid


def mask_wall(df):
    df = df.copy()
    is_wall = ((df["x"] < X_CORNER) & (df["y"] >= Y_CORNER))
    if "rho" in df.columns:
        is_wall = is_wall | (df["rho"] > WALL_RHO)
    for col in ["rho", "u", "v", "P", "e", "mf"]:
        if col in df.columns:
            df.loc[is_wall, col] = np.nan
    return df


def add_wall(ax):
    ax.add_patch(plt.Rectangle((0, Y_CORNER), X_CORNER, Y_MAX - Y_CORNER,
                                facecolor="#888888", edgecolor="black",
                                linewidth=1.5, zorder=5))


def setup_axes(ax):
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def save_pressure_frame(f, t, pics_path, idx):
    """Сохраняет один кадр давления в отдельный файл pressure_tN.png."""
    df = read_csv(f)
    if df is None: return
    df_t = mask_wall(df[df["t"] == df["t"].max()])
    xs, ys, grid = to_grid(df_t, "P")

    fig, ax = plt.subplots(figsize=(ASPECT * 3.0, 3.0))
    extent = [xs.min(), xs.max(), ys.min(), ys.max()]
    im = ax.imshow(grid, origin="lower", extent=extent,
                   aspect="equal", cmap="hot_r",
                   vmin=0, vmax=P_DET * 1.5)
    plt.colorbar(im, ax=ax, label="P", fraction=0.025, pad=0.02)
    setup_axes(ax)
    add_wall(ax)
    ax.set_title(f"t = {t:.3f}")

    plt.tight_layout()
    out = os.path.join(pics_path, f"pressure_t{idx}.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print("Сохранено:", out)


def plot_pressure_snapshots(csv_folder, pics_path):
    step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
    if not step_files:
        print("Нет step-файлов"); return []

    records = []
    for f in step_files:
        df = read_csv(f)
        if df is None or df.empty: continue
        t = df["t"].max()
        df_t = df[df["t"] == t]
        front = df_t[df_t["P"] > P_THRESH]
        xf = front["x"].max() if not front.empty else 0.0
        records.append((t, xf, f))

    if not records: return []

    # Три кадра: (1) фронт в канале, (2) момент выхода, (3) после дифракции
    records.sort(key=lambda r: r[0])

    # t1: где x_front ближе всего к середине канала
    r1 = min(records, key=lambda r: abs(r[1] - X_CORNER * 0.5))
    # t2: где x_front ближе всего к концу канала
    r2 = min(records, key=lambda r: abs(r[1] - X_CORNER))
    # t3: последний момент (дифракция)
    r3 = records[-1]

    chosen = [r1, r2, r3]
    for i, (t, xf, f) in enumerate(chosen, 1):
        save_pressure_frame(f, t, pics_path, i)

    return chosen


def plot_mass_fraction_map(csv_folder, pics_path):
    final = os.path.join(csv_folder, "Final.csv")
    if not os.path.exists(final):
        step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
        if not step_files: print("Нет финального файла"); return None
        final = step_files[-1]

    df = read_csv(final)
    if df is None or "mf" not in df.columns: return None

    t = df["t"].max()
    df_t = mask_wall(df[df["t"] == t])
    xs, ys, grid = to_grid(df_t, "mf")

    fig, ax = plt.subplots(figsize=(ASPECT * 3.0, 3.0))
    extent = [xs.min(), xs.max(), ys.min(), ys.max()]
    im = ax.imshow(grid, origin="lower", extent=extent,
                   aspect="equal", cmap="YlOrBr", vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="W", fraction=0.025, pad=0.02)
    setup_axes(ax)
    add_wall(ax)
    ax.set_title(f"t = {t:.3f}")

    plt.tight_layout()
    out = os.path.join(pics_path, "W_map.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print("Сохранено:", out)
    return df_t


def plot_arrival_time(csv_folder, pics_path):
    step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
    if not step_files:
        print("Нет step-файлов для t_arrival"); return None

    arrival = {}
    for f in step_files:
        df = read_csv(f)
        if df is None or df.empty: continue
        t = df["t"].max()
        front = df[(df["t"] == t) & (df["P"] > P_THRESH)]
        for _, row in front.iterrows():
            key = (row["x"], row["y"])
            if key not in arrival:
                arrival[key] = t

    if not arrival:
        print("Фронт не обнаружен"); return None

    xs_arr = np.array([k[0] for k in arrival])
    ys_arr = np.array([k[1] for k in arrival])
    ta_arr = np.array([arrival[k] for k in arrival])

    xs_u = np.sort(np.unique(xs_arr))
    ys_u = np.sort(np.unique(ys_arr))
    grid_ta = np.full((len(ys_u), len(xs_u)), np.nan)

    xi_idx = {v: i for i, v in enumerate(xs_u)}
    yi_idx = {v: i for i, v in enumerate(ys_u)}
    for k in range(len(xs_arr)):
        xi = xi_idx.get(xs_arr[k])
        yi = yi_idx.get(ys_arr[k])
        if xi is not None and yi is not None:
            grid_ta[yi, xi] = ta_arr[k]

    XX, YY = np.meshgrid(xs_u, ys_u)
    grid_ta[(XX < X_CORNER) & (YY < Y_CORNER) & np.isnan(grid_ta)] = 0.0

    t_max = np.nanmax(grid_ta)
    t_at_corner = X_CORNER / D_CJ
    all_levels = np.linspace(t_at_corner, t_max, 3)

    fig, ax = plt.subplots(figsize=(ASPECT * 3.0, 3.0))
    extent = [xs_u.min(), xs_u.max(), ys_u.min(), ys_u.max()]

    im = ax.imshow(grid_ta, origin="lower", extent=extent,
                   aspect="equal", cmap="viridis",
                   vmin=0, vmax=t_max)
    plt.colorbar(im, ax=ax, label="t", fraction=0.025, pad=0.02)

    cs = ax.contour(xs_u, ys_u, grid_ta,
                    levels=all_levels,
                    colors="white", linewidths=0.9, alpha=0.85)
    ax.clabel(cs, fmt="%.2f", fontsize=7, colors="white", inline=True)

    setup_axes(ax)
    add_wall(ax)

    plt.tight_layout()
    out = os.path.join(pics_path, "arrival_time.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print("Сохранено:", out)
    return (xs_u, ys_u, grid_ta)


def write_values_mader3(csv_folder, pics_path):
    """4c: контрольные числа для corner turning."""
    step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
    if not step_files: return

    # Скорость фронта в канале (y в нижней полосе, x < X_CORNER)
    # берём точки, где фронт ещё в канале
    ts_in_ch, xfs_in_ch = [], []
    for f in step_files:
        df = read_csv(f)
        if df is None or df.empty: continue
        t = df["t"].max()
        df_t = df[df["t"] == t]
        # центральный ряд канала
        ys = np.sort(df_t[df_t["y"] < Y_CORNER]["y"].unique())
        if len(ys) == 0: continue
        y_mid = ys[len(ys) // 2]
        row = df_t[np.isclose(df_t["y"], y_mid)]
        front = row[(row["P"] > P_THRESH) & (row["x"] < X_CORNER)]
        if front.empty: continue
        xf = front["x"].max()
        if xf < X_CORNER - 0.1:   # фронт ещё не у выхода
            ts_in_ch.append(t)
            xfs_in_ch.append(xf)


def main():
    if len(sys.argv) < 2:
        print("Использование: python3 PlotMader3.py <config_name>"); return

    config_name = sys.argv[1]
    csv_folder  = os.path.join("output", config_name, "CSV")
    pics_path   = os.path.join("results", "04_mader_corner")

    if not os.path.exists(csv_folder):
        print("Папка CSV не найдена:", csv_folder); return

    os.makedirs(pics_path, exist_ok=True)

    plot_pressure_snapshots(csv_folder, pics_path)   # pressure_t1/t2/t3.png
    plot_mass_fraction_map(csv_folder, pics_path)    # W_map.png
    plot_arrival_time(csv_folder, pics_path)         # arrival_time.png
    write_values_mader3(csv_folder, pics_path)       # values.txt




if __name__ == "__main__":
    main()