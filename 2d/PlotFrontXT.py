# PlotFrontXT.py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys, os, glob

plt.rcParams.update({"font.size": 12, "figure.dpi": 120})

GAMMA = 3.0
RHO_0 = 1.84
D_CJ  = 0.88
X_DET = 0.03
P_CJ  = RHO_0 * D_CJ**2 / (GAMMA + 1)
P_THRESH = 0.05 * P_CJ

# Цвета для разных сеток
GRID_COLORS = ["#1f77b4", "#2ca02c", "#9467bd", "#8c564b", "#e377c2"]


def collect_front(csv_folder):
    """Собирает массивы ts, xfs по всем step-файлам расчёта."""
    step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
    if not step_files:
        return None, None

    ts, xfs = [], []
    for f in step_files:
        df = pd.read_csv(f)
        if df.empty: continue
        t = df["t"].max()
        df_t = df[df["t"] == t]
        ys = np.sort(df_t["y"].unique())
        y_mid = ys[len(ys) // 2]
        row = df_t[np.isclose(df_t["y"], y_mid)]
        front = row[row["P"] > P_THRESH]
        if front.empty: continue
        ts.append(t)
        xfs.append(front["x"].max())

    if not ts:
        return None, None
    return np.array(ts), np.array(xfs)


def fit_slope(ts, xfs):
    """Линейная регрессия по последним 70% точек."""
    n = len(ts)
    start = max(1, int(0.3 * n))
    slope, _ = np.polyfit(ts[start:], xfs[start:], 1)
    err = 100 * abs(slope - D_CJ) / D_CJ
    return slope, err


def detect_grid_label(csv_folder):
    """Определяет N_x из первого CSV-файла для подписи легенды."""
    step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
    if not step_files:
        step_files = [os.path.join(csv_folder, "Final.csv")]
    try:
        df = pd.read_csv(step_files[0])
        t0 = df["t"].max()
        df = df[df["t"] == t0]
        ys = np.sort(df["y"].unique())
        y_mid = ys[len(ys) // 2]
        row = df[np.isclose(df["y"], y_mid)]
        return len(row["x"].unique())
    except Exception:
        return None


def plot_front_xt_single(ts, xfs, out_path):
    """Один расчёт — одна кривая + аналитика."""
    slope, err = fit_slope(ts, xfs)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ts, xfs, "-", color="#1f77b4", linewidth=2.0,
            label=f"Численное решение (D = {slope:.4f}, err = {err:.2f}%)")
    t_line = np.linspace(ts[0], ts[-1], 100)
    ax.plot(t_line, X_DET + D_CJ * t_line, "--", color="#d62728",
            linewidth=1.5,
            label=f"Аналитическое решение (D = $D_{{CJ}}$ = {D_CJ:.4f})")

    ax.set_xlabel("t"); ax.set_ylabel("x")
    ax.set_title("Положение фронта детонации")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print("Сохранено:", out_path)
    return slope, err


def plot_front_xt_multi(runs, out_path):
    """Несколько расчётов — несколько кривых + одна аналитика.
    runs: список (label, ts, xfs, slope, err)."""
    fig, ax = plt.subplots(figsize=(9, 5.5))

    t_all_min = min(r[1].min() for r in runs)
    t_all_max = max(r[1].max() for r in runs)

    for i, (label, ts, xfs, slope, err) in enumerate(runs):
        color = GRID_COLORS[i % len(GRID_COLORS)]
        ax.plot(ts, xfs, "-", color=color, linewidth=2.0,
                label=f"{label} (D = {slope:.4f}, err = {err:.2f}%)")

    t_line = np.linspace(t_all_min, t_all_max, 100)
    ax.plot(t_line, X_DET + D_CJ * t_line, "--", color="#d62728", linewidth=1.5,
            label=f"Аналитическое решение (D = $D_{{CJ}}$ = {D_CJ:.4f})")

    ax.set_xlabel("t"); ax.set_ylabel("x")
    ax.set_title("Положение фронта детонации")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print("Сохранено:", out_path)


def plot_profiles(csv_folder, out_path):
    final = os.path.join(csv_folder, "Final.csv")
    if not os.path.exists(final):
        step_files = sorted(glob.glob(os.path.join(csv_folder, "*_step.csv")))
        if not step_files: return None
        final = step_files[-1]

    df = pd.read_csv(final)
    if df.empty: return None

    t_max = df["t"].max()
    df = df[df["t"] == t_max]

    ys = np.sort(df["y"].unique())
    y_mid = ys[len(ys) // 2]
    row = df[np.isclose(df["y"], y_mid)].sort_values("x")

    x_arr   = row["x"].values
    P_arr   = row["P"].values
    rho_arr = row["rho"].values
    has_mf  = "mf" in row.columns
    W_arr   = row["mf"].values if has_mf else None

    x_front_th = X_DET + D_CJ * t_max
    rho_CJ = RHO_0 * (GAMMA + 1) / GAMMA

    fig, axs = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    axs[0].plot(x_arr, P_arr, color="#2ca02c", linewidth=2.0, label="Численное решение")
    axs[0].set_ylabel(r"$P$"); axs[0].legend(loc="best", fontsize=10)
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(x_arr, rho_arr, color="#1f77b4", linewidth=2.0, label="Численное решение")
    axs[1].set_ylabel(r"$\rho$"); axs[1].legend(loc="best", fontsize=10)
    axs[1].grid(True, alpha=0.3)

    if has_mf:
        axs[2].plot(x_arr, W_arr, color="#8c564b", linewidth=2.0, label="Численное решение")
        axs[2].set_ylabel(r"$W$"); axs[2].legend(loc="best", fontsize=10)
        axs[2].grid(True, alpha=0.3); axs[2].set_ylim(-0.05, 1.1)
    axs[2].set_xlabel("x")

    fig.suptitle(f"Профили P, ρ, W  (t = {t_max:.4f})", fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print("Сохранено:", out_path)

    return float(np.max(P_arr))


def main():
    if len(sys.argv) < 2:
        print("Использование:")
        print("  Один расчёт:          python3 PlotFrontXT.py <config_name>")
        print("  Несколько сеток:      python3 PlotFrontXT.py <cfg_N100> <cfg_N200> <cfg_N400>")
        return

    config_names = sys.argv[1:]
    results_dir  = os.path.join("results", "04_mader_det")
    os.makedirs(results_dir, exist_ok=True)

    # --- Собираем данные по всем расчётам ---
    runs = []
    for cfg in config_names:
        csv_folder = os.path.join("output", cfg, "CSV")
        if not os.path.isdir(csv_folder):
            print(f"Пропускаю {cfg}: нет папки {csv_folder}")
            continue
        ts, xfs = collect_front(csv_folder)
        if ts is None:
            print(f"Пропускаю {cfg}: фронт не обнаружен")
            continue
        slope, err = fit_slope(ts, xfs)
        N = detect_grid_label(csv_folder)
        label = f"N = {N}" if N else cfg
        runs.append((label, ts, xfs, slope, err, cfg))

    if not runs:
        print("Нет валидных расчётов"); return

    # --- front_xt.png ---
    if len(runs) == 1:
        # один расчёт: сохраняем и в pics/, и в results/
        cfg = runs[0][5]
        pics_path = os.path.join("output", cfg, "pics")
        os.makedirs(pics_path, exist_ok=True)
        ts, xfs = runs[0][1], runs[0][2]
        plot_front_xt_single(ts, xfs, os.path.join(pics_path, "front_xt.png"))
        plot_front_xt_single(ts, xfs, os.path.join(results_dir, "front_xt.png"))
    else:
        # несколько расчётов: один общий график в results/
        # (без поля runs[i][5] — срезаем последнее для сигнатуры)
        runs_for_plot = [(r[0], r[1], r[2], r[3], r[4]) for r in runs]
        plot_front_xt_multi(runs_for_plot, os.path.join(results_dir, "front_xt.png"))

    # --- profiles_P_rho_W.png: берём ПЕРВЫЙ расчёт (базовый) ---
    base_run = next((r for r in runs if "200" in r[0]), runs[0])
    base_cfg = base_run[5]
    csv_folder = os.path.join("output", base_cfg, "CSV")
    pics_path  = os.path.join("output", base_cfg, "pics")
    os.makedirs(pics_path, exist_ok=True)
    P_max = plot_profiles(csv_folder, os.path.join(pics_path, "profiles_P_rho_W.png"))
    plot_profiles(csv_folder, os.path.join(results_dir, "profiles_P_rho_W.png"))

    # --- values.txt ---
    err_P = 100 * abs(P_max - P_CJ) / P_CJ if P_max is not None else 0.0
    vals_path = os.path.join(results_dir, "values.txt")
    with open(vals_path, "w") as f:
        f.write("Мейдер 2: Одномерная детонация CJ\n\n")
        for (label, _, _, slope, err, _) in runs:
            f.write(f"{label}: D_computed = {slope:.4f}, "
                    f"D_CJ = {D_CJ:.4f}, err = {err:.2f}%\n")
        if P_max is not None:
            f.write(f"P_max = {P_max:.4f}, P_CJ = {P_CJ:.4f}, err = {err_P:.2f}%\n")
    print("Сохранено:", vals_path)


if __name__ == "__main__":
    main()