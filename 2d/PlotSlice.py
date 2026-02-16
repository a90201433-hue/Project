import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os


plt.style.use("seaborn-v0_8-whitegrid")

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.dpi": 120,
})


def read_data(filename):
    try:
        return pd.read_csv(filename)
    except Exception:
        print(f"Ошибка чтения файла: {filename}")
        return None



def main():

    if len(sys.argv) < 2:
        print("Использование:")
        print("python3 PlotSlice.py <config_name> [y_value]")
        return

    config_name = sys.argv[1]


    if len(sys.argv) >= 3:
        fixed_y = float(sys.argv[2])
    else:
        fixed_y = None


    base_path = os.path.join("output", config_name)
    csv_path = os.path.join(base_path, "CSV", "Final.csv")
    pics_path = os.path.join(base_path, "pics")

    if not os.path.exists(csv_path):
        print("Файл Final.csv не найден:", csv_path)
        return

    os.makedirs(pics_path, exist_ok=True)

    data = read_data(csv_path)
    if data is None:
        return

    t_max = data["t"].max()
    data = data[data["t"] == t_max]


    if fixed_y is None:
        y_vals = np.sort(data["y"].unique())
        fixed_y = y_vals[len(y_vals)//2]

    slice_df = data[np.isclose(data["y"], fixed_y)]
    slice_df = slice_df.sort_values("x")

    if slice_df.empty:
        print(f"Нет данных для y = {fixed_y}")
        return


    # --- Приятная палитра ---
    colors = {
        "rho": "#1f77b4",   # мягкий синий
        "u":   "#d62728",   # спокойный красный
        "P":   "#2ca02c",   # зелёный
        "e":   "#9467bd"    # фиолетовый
    }

    fig, axs = plt.subplots(2, 2, figsize=(9, 7), sharex=True)

    lw = 2.2

    axs[0, 0].plot(slice_df["x"], slice_df["rho"],
                   color=colors["rho"], linewidth=lw)
    axs[0, 0].set_ylabel(r'$\rho$')

    axs[0, 1].plot(slice_df["x"], slice_df["u"],
                   color=colors["u"], linewidth=lw)
    axs[0, 1].set_ylabel(r'$u$')

    axs[1, 0].plot(slice_df["x"], slice_df["P"],
                   color=colors["P"], linewidth=lw)
    axs[1, 0].set_ylabel(r'$P$')

    axs[1, 1].plot(slice_df["x"], slice_df["e"],
                   color=colors["e"], linewidth=lw)
    axs[1, 1].set_ylabel(r'$e$')

    for ax in axs.flatten():
        ax.set_xlabel(r'$x$')
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(
        rf"{config_name}   |   $y={fixed_y:.4f}$   |   $t={t_max:.4f}$",
        fontsize=14
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    output_image = os.path.join(pics_path, "Final.png")
    plt.savefig(output_image, dpi=300)
    plt.close()

    print("Сохранено:", output_image)



if __name__ == "__main__":
    main()
