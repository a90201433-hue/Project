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
        print("python3 PlotSlice.py <config_name> [y_value] [csv_name]")
        return

    config_name = sys.argv[1]

    fixed_y = None
    csv_name = "Final.csv"   # по умолчанию

    # --- разбор аргументов ---
    if len(sys.argv) >= 3:
        try:
            fixed_y = float(sys.argv[2])
        except ValueError:
            csv_name = sys.argv[2]

    if len(sys.argv) >= 4:
        csv_name = sys.argv[3]

    base_path = os.path.join("output", config_name)
    csv_path = os.path.join(base_path, "CSV", csv_name)
    pics_path = os.path.join(base_path, "pics")

    if not os.path.exists(csv_path):
        print("CSV файл не найден:", csv_path)
        return

    os.makedirs(pics_path, exist_ok=True)

    data = read_data(csv_path)
    if data is None:
        return

    # берём последний момент времени
    t_max = data["t"].max()
    data = data[data["t"] == t_max]

    # если y не задан — берём центральный
    if fixed_y is None:
        y_vals = np.sort(data["y"].unique())
        fixed_y = y_vals[len(y_vals)//2]

    slice_df = data[np.isclose(data["y"], fixed_y)]
    slice_df = slice_df.sort_values("x")

    if slice_df.empty:
        print(f"Нет данных для y = {fixed_y}")
        return

    # --- цвета ---
    colors = {
        "rho": "#1f77b4",
        "u":   "#d62728",
        "P":   "#2ca02c",
        "e":   "#9467bd"
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
        rf"{config_name} | $y={fixed_y:.4f}$ | $t={t_max:.4f}$",
        fontsize=14
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # имя картинки = имя csv без расширения
    image_name = os.path.splitext(csv_name)[0] + ".png"
    output_image = os.path.join(pics_path, image_name)

    plt.savefig(output_image, dpi=300)
    plt.close()

    print("Сохранено:", output_image)


if __name__ == "__main__":
    main()