import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.optimize import brentq


plt.style.use("seaborn-v0_8-whitegrid")

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.dpi": 120,
})


def exact_wave_positions_sod(gamma=1.4, t=0.2, x_gap=0.5):
    rho_L, u_L, P_L = 1.0, 0.0, 1.0
    rho_R, u_R, P_R = 0.125, 0.0, 0.1
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)

    def f_K(P, rho_K, P_K, c_K):
        if P > P_K:
            A_K = 2.0 / ((gamma + 1) * rho_K)
            B_K = (gamma - 1) / (gamma + 1) * P_K
            return (P - P_K) * np.sqrt(A_K / (P + B_K))
        return (2 * c_K / (gamma - 1)) * ((P / P_K) ** ((gamma - 1) / (2 * gamma)) - 1)

    P_star = brentq(lambda P: f_K(P, rho_L, P_L, c_L)
                            + f_K(P, rho_R, P_R, c_R) + (u_R - u_L),
                    1e-12, max(P_L, P_R) * 10)
    u_star = 0.5 * (u_L + u_R) + 0.5 * (f_K(P_star, rho_R, P_R, c_R)
                                       - f_K(P_star, rho_L, P_L, c_L))

    S_HL = u_L - c_L
    c_sL = c_L * (P_star / P_L) ** ((gamma - 1) / (2 * gamma))
    S_TL = u_star - c_sL
    S_CR = u_star
    S_R = u_R + c_R * np.sqrt((gamma + 1) / (2 * gamma) * P_star / P_R
                              + (gamma - 1) / (2 * gamma))

    return {
        "VR_head": x_gap + S_HL * t,
        "VR_tail": x_gap + S_TL * t,
        "KR":      x_gap + S_CR * t,
        "UV":      x_gap + S_R  * t,
    }


def annotate_waves(ax, positions):
    y_min, y_max = ax.get_ylim()
    y_label = y_max - 0.07 * (y_max - y_min)

    line_color = "#555555"
    text_color = "#333333"
    text_bbox  = dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="none", alpha=0.75)

    ax.axvline(positions["VR_head"], color=line_color, linestyle="--",
               linewidth=0.9, alpha=0.8, zorder=1)
    ax.axvline(positions["VR_tail"], color=line_color, linestyle="--",
               linewidth=0.9, alpha=0.8, zorder=1)
    ax.text(0.5 * (positions["VR_head"] + positions["VR_tail"]), y_label,
            "ВР", ha="center", va="top", fontsize=9, color=text_color,
            bbox=text_bbox)

    ax.axvline(positions["KR"], color=line_color, linestyle="--",
               linewidth=0.9, alpha=0.8, zorder=1)
    ax.text(positions["KR"], y_label, "КР", ha="center", va="top",
            fontsize=9, color=text_color, bbox=text_bbox)

    ax.axvline(positions["UV"], color=line_color, linestyle="--",
               linewidth=0.9, alpha=0.8, zorder=1)
    ax.text(positions["UV"], y_label, "УВ", ha="center", va="top",
            fontsize=9, color=text_color, bbox=text_bbox)



def read_data(filename):
    try:
        return pd.read_csv(filename)
    except Exception:
        print(f"Ошибка чтения файла: {filename}")
        return None


def save_sod_single_plots(slice_df, coord, t_max, results_dir, colors, lw):
    """4a: три отдельных графика rho, P, u в results/04_mader_sod/."""
    os.makedirs(results_dir, exist_ok=True)
    wave_pos = exact_wave_positions_sod(gamma=1.4, t=t_max, x_gap=0.5)
    x_arr = slice_df[coord].values

    plots = [
        ("sod_rho.png", "rho", slice_df["rho"].values, r'$\rho$'),
        ("sod_p.png",   "P",   slice_df["P"].values,   r'$P$'),
        ("sod_u.png",   "vel", slice_df["u"].values,   r'$u$'),
    ]

    for fname, color_key, y_num, ylabel in plots:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(x_arr, y_num, color=colors[color_key], linewidth=lw,
                label="Численное решение")
        ax.set_xlabel(f'${coord}$')
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(f"Мейдер 1: тест Сода  $t={t_max:.3f}$")

        annotate_waves(ax, wave_pos)
        ax.legend(loc="best", fontsize=10)

        plt.tight_layout()
        out = os.path.join(results_dir, fname)
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print("Сохранено:", out)


def write_values_mader1(slice_df, coord_col, pics_path):
    vals_path = os.path.join(pics_path, "values.txt")
    ref = {
        0.60: {"rho": 0.4263, "p": 0.3031, "u": 0.9274},
        0.75: {"rho": 0.2656, "p": 0.3031, "u": 0.9274},
    }
    with open(vals_path, "w") as f:
        f.write("Мейдер 1: Ударная труба Сода \n\n")
        for x_target, r in ref.items():
            idx = (slice_df[coord_col] - x_target).abs().idxmin()
            row = slice_df.loc[idx]
            rho, p, u = row["rho"], row["P"], row["u"]
            err_rho = 100 * abs(rho - r["rho"]) / r["rho"]
            err_p   = 100 * abs(p   - r["p"])   / r["p"]
            err_u   = 100 * abs(u   - r["u"])   / r["u"]
            f.write(f"x = {x_target:.2f}:\n")
            f.write(f"rho = {rho:.4f}, ref_rho = {r['rho']:.4f}, err = {err_rho:.2f}%\n")
            f.write(f"p = {p:.4f}, ref_p = {r['p']:.4f}, err = {err_p:.2f}%\n")
            f.write(f"u = {u:.4f}, ref_u = {r['u']:.4f}, err = {err_u:.2f}%\n\n")
    print("Сохранено:", vals_path)


def write_values_mader2(slice_df, coord_col, pics_path, t_max):
    vals_path = os.path.join(pics_path, "values.txt")

    gamma, rho_0, D_CJ = 3.0, 1.84, 0.88
    P_CJ = rho_0 * D_CJ**2 / (gamma + 1)   # 0.3562
    x_det = 0.03
    P_threshold = 0.05 * P_CJ

    front = slice_df[slice_df["P"] > P_threshold]
    x_front = front[coord_col].max() if not front.empty else x_det
    D_computed = (x_front - x_det) / t_max if t_max > 0 else 0.0

    P_max = slice_df["P"].max()

    err_D = 100 * abs(D_computed - D_CJ) / D_CJ
    err_P = 100 * abs(P_max - P_CJ) / P_CJ

    with open(vals_path, "w") as f:
        f.write("Мейдер 2: Одномерная детонация CJ\n\n")
        f.write(f"D_computed = {D_computed:.4f}, D_CJ = {D_CJ:.4f}, err = {err_D:.2f}%\n")
        f.write(f"P_max = {P_max:.4f}, P_CJ = {P_CJ:.4f}, err = {err_P:.2f}%\n")
    print("Сохранено:", vals_path)

def main():
    if len(sys.argv) < 2:
        print("Использование: python3 PlotSlice.py <config_name> [axis] [value] [csv_name]")
        return

    config_name = sys.argv[1]
    axis        = "y"
    fixed_value = None
    csv_name    = "Final.csv"
    args        = sys.argv[2:]

    if len(args) >= 1 and args[0] in ["x", "y"]:
        axis = args[0]; args = args[1:]
    if len(args) >= 1:
        try: fixed_value = float(args[0]); args = args[1:]
        except ValueError: pass
    if len(args) >= 1:
        csv_name = args[0]

    base_path = os.path.join("output", config_name)
    csv_path  = os.path.join(base_path, "CSV", csv_name)
    pics_path = os.path.join(base_path, "pics")

    if not os.path.exists(csv_path):
        print("CSV файл не найден:", csv_path); return

    os.makedirs(pics_path, exist_ok=True)

    data  = read_data(csv_path)
    if data is None: return

    t_max = data["t"].max()
    data  = data[data["t"] == t_max]

    if fixed_value is None:
        unique_vals = np.sort(data[axis].unique())
        fixed_value = unique_vals[len(unique_vals) // 2]

    slice_df = data[np.isclose(data[axis], fixed_value)]
    if slice_df.empty:
        print(f"Нет данных для {axis} = {fixed_value}"); return

    coord              = "x"   if axis == "y" else "y"
    velocity_component = "u"   if axis == "y" else "v"
    velocity_label     = r"$u$" if axis == "y" else r"$v$"

    slice_df = slice_df.sort_values(coord)
    x_arr    = slice_df[coord].values
    x_arr_th = np.linspace(x_arr.min(), x_arr.max(), 2000)

    has_mf = "mf" in slice_df.columns

    colors = {
        "rho": "#1f77b4", "vel": "#d62728",
        "P":   "#2ca02c", "e":   "#9467bd",
        "mf":  "#8c564b", "ana": "#ff7f0e",
        "ref": "#555555",
    }

    mader2_refs   = None
    mader2_xfront = None

    if config_name == "Mader2":
        gamma  = 3.0; rho_0 = 1.84; D_CJ = 0.88; x_det = 0.03
        P_CJ   = rho_0 * D_CJ**2 / (gamma + 1)
        rho_CJ = rho_0 * (gamma + 1) / gamma
        u_CJ   = D_CJ / (gamma + 1)
        e_CJ   = P_CJ / ((gamma - 1) * rho_CJ)
        mader2_xfront = x_det + D_CJ * t_max
        mader2_refs   = {
            "rho": (rho_CJ, r"$\rho_{CJ} = %.4f$" % rho_CJ),
            "u":   (u_CJ,   r"$u_{CJ} = %.4f$"   % u_CJ),
            "P":   (P_CJ,   r"$P_{CJ} = %.4f$"   % P_CJ),
            "e":   (e_CJ,   r"$e_{CJ} = %.4f$"   % e_CJ),
        }


    is_mader2 = (config_name == "Mader2") and has_mf

    if is_mader2:
        fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    else:
        fig, axs = plt.subplots(2, 2, figsize=(9, 7), sharex=True)

    lw     = 2.2

    def draw(ax, key, color_key, y_num, y_ana, ylabel, ref_key=None):
        ax.plot(x_arr, y_num, color=colors[color_key], linewidth=lw, label="Численное решение")

        if mader2_refs is not None and ref_key is not None:
            ax.axvline(mader2_xfront, color=colors["ref"], linewidth=1.3,
                    linestyle=":", label=f"$x_{{D_{{CJ}}}}={mader2_xfront:.3f}$")
            ax.legend(fontsize=8, loc="best")

        ax.set_ylabel(ylabel)
        ax.set_xlabel(f'${coord}$')
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    if is_mader2:
        draw(axs[0,0], "rho", "rho", slice_df["rho"].values,               None, r'$\rho$',   "rho")
        draw(axs[0,1], "u",   "vel", slice_df[velocity_component].values,  None, velocity_label, "u")
        draw(axs[0,2], "P",   "P",   slice_df["P"].values,                 None, r'$P$',      "P")
        draw(axs[1,0], "e",   "e",   slice_df["e"].values,                 None, r'$e$',      "e")

        
        ax_mf = axs[1, 1]
        ax_mf.plot(x_arr, slice_df["mf"].values,
                   color=colors["mf"], linewidth=lw, label="Численное решение")
        mf_theory = np.where(x_arr_th < mader2_xfront, 0.0, 1.0)
        ax_mf.axvline(mader2_xfront, color=colors["ref"], linewidth=1.3,
                      linestyle=":", label=f"$x_{{D_{{CJ}}}}={mader2_xfront:.3f}$")
        ax_mf.set_ylabel(r'$W$')
        ax_mf.set_xlabel(f'${coord}$')
        ax_mf.legend(fontsize=8, loc="best")
        ax_mf.grid(True, alpha=0.3)
        ax_mf.spines['top'].set_visible(False)
        ax_mf.spines['right'].set_visible(False)

        axs[1, 2].set_visible(False)

    else:
        draw(axs[0,0], "rho", "rho", slice_df["rho"].values,               None, r'$\rho$')
        draw(axs[0,1], "u",   "vel", slice_df[velocity_component].values,  None, velocity_label)
        draw(axs[1,0], "P",   "P",   slice_df["P"].values,                 None, r'$P$')
        draw(axs[1,1], "e",   "e",   slice_df["e"].values,                 None, r'$e$')


    if config_name == "Mader1":
        wave_pos = exact_wave_positions_sod(gamma=1.4, t=t_max, x_gap=0.5)
        for ax in [axs[0,0], axs[0,1], axs[1,0], axs[1,1]]:
            annotate_waves(ax, wave_pos)

    fig.suptitle(
        rf"{config_name} | {axis}={fixed_value:.4f} | $t={t_max:.4f}$",
        fontsize=14
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    image_name   = os.path.splitext(csv_name)[0] + f"_{axis}.png"
    output_image = os.path.join(pics_path, image_name)
    plt.savefig(output_image, dpi=300)
    plt.close()

    print("Сохранено:", output_image)

    if config_name == "Mader1":
        write_values_mader1(slice_df, coord, pics_path)
        results_dir = os.path.join("results", "04_mader_sod")
        save_sod_single_plots(slice_df, coord, t_max, results_dir, colors, lw)
        write_values_mader1(slice_df, coord, results_dir)



if __name__ == "__main__":
    main()