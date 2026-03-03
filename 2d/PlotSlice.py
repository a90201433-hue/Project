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
        print("python3 PlotSlice.py <config_name> [axis] [value] [csv_name]")
        return

    config_name = sys.argv[1]

    axis = "y"            # по умолчанию
    fixed_value = None
    csv_name = "Final.csv"

    args = sys.argv[2:]

    # --- axis ---
    if len(args) >= 1 and args[0] in ["x", "y"]:
        axis = args[0]
        args = args[1:]

    # --- значение ---
    if len(args) >= 1:
        try:
            fixed_value = float(args[0])
            args = args[1:]
        except ValueError:
            pass

    # --- имя csv ---
    if len(args) >= 1:
        csv_name = args[0]

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

    # --- последний момент времени ---
    t_max = data["t"].max()
    data = data[data["t"] == t_max]

    # --- если значение не задано → берём центральное ---
    if fixed_value is None:
        unique_vals = np.sort(data[axis].unique())
        fixed_value = unique_vals[len(unique_vals)//2]

    # --- формируем сечение ---
    slice_df = data[np.isclose(data[axis], fixed_value)]

    if slice_df.empty:
        print(f"Нет данных для {axis} = {fixed_value}")
        return

    # --- определяем направление построения ---
    if axis == "y":
        coord = "x"
        velocity_component = "u"
        velocity_label = r"$u$"
    else:
        coord = "y"
        velocity_component = "v"
        velocity_label = r"$v$"

    slice_df = slice_df.sort_values(coord)

    # --- цвета ---
    colors = {
        "rho": "#1f77b4",
        "vel": "#d62728",
        "P":   "#2ca02c",
        "e":   "#9467bd"
    }

    fig, axs = plt.subplots(2, 2, figsize=(9, 7), sharex=True)

    lw = 2.2

    # rho
    axs[0, 0].plot(slice_df[coord], slice_df["rho"],
                   color=colors["rho"], linewidth=lw)
    axs[0, 0].set_ylabel(r'$\rho$')

    # velocity
    axs[0, 1].plot(slice_df[coord], slice_df[velocity_component],
                   color=colors["vel"], linewidth=lw)
    axs[0, 1].set_ylabel(velocity_label)

    # pressure
    axs[1, 0].plot(slice_df[coord], slice_df["P"],
                   color=colors["P"], linewidth=lw)
    axs[1, 0].set_ylabel(r'$P$')

    # internal energy
    axs[1, 1].plot(slice_df[coord], slice_df["e"],
                   color=colors["e"], linewidth=lw)
    axs[1, 1].set_ylabel(r'$e$')

    for ax in axs.flatten():
        ax.set_xlabel(f'${coord}$')
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(
        rf"{config_name} | {axis}={fixed_value:.4f} | $t={t_max:.4f}$",
        fontsize=14
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    image_name = os.path.splitext(csv_name)[0] + f"_{axis}.png"
    output_image = os.path.join(pics_path, image_name)

    plt.savefig(output_image, dpi=300)
    plt.close()

    print("Сохранено:", output_image)


if __name__ == "__main__":
    main()