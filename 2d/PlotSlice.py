import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

plt.style.use('seaborn-v0_8')


# --- ЧТЕНИЕ ДАННЫХ ---

def read_data(filename):
    try:
        return pd.read_csv(filename)
    except Exception:
        print(f"Ошибка чтения файла: {filename}")
        return None


# --- ОСНОВНАЯ ЛОГИКА ---

def main():

    # --- Аргументы командной строки ---
    if len(sys.argv) < 2:
        print("Использование:")
        print("python3 plot_slice.py <csv_file> [y_value]")
        return

    filename = sys.argv[1]

    if len(sys.argv) >= 3:
        fixed_y = float(sys.argv[2])
    else:
        fixed_y = None

    data = read_data(filename)
    if data is None:
        return

    # Берём последнее время
    t_max = data["t"].max()
    data = data[data["t"] == t_max]

    # Выбор среза по y
    if fixed_y is None:
        y_vals = np.sort(data["y"].unique())
        fixed_y = y_vals[len(y_vals)//2]

    slice_df = data[np.isclose(data["y"], fixed_y)]
    slice_df = slice_df.sort_values("x")

    if slice_df.empty:
        print(f"Нет данных для y = {fixed_y}")
        return

    # --- Создание фигуры ---
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

    axs[0, 0].plot(slice_df["x"], slice_df["rho"], color='blue')
    axs[0, 0].set_ylabel(r'$\rho$', fontsize=14)

    axs[0, 1].plot(slice_df["x"], slice_df["u"], color='red')
    axs[0, 1].set_ylabel('u', fontsize=14)

    axs[1, 0].plot(slice_df["x"], slice_df["P"], color='green')
    axs[1, 0].set_ylabel('P', fontsize=14)

    axs[1, 1].plot(slice_df["x"], slice_df["e"], color='purple')
    axs[1, 1].set_ylabel(r'$e$', fontsize=14)

    for ax in axs.flatten():
        ax.grid(True, alpha=0.5)
        ax.set_xlabel('x', fontsize=14)
        ax.tick_params(labelsize=12)

    fig.suptitle(f"Slice at y = {fixed_y:.6f}, t = {t_max:.6f}", fontsize=16)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # --- Автоматическое имя картинки ---
    base = os.path.splitext(os.path.basename(filename))[0]
    output_image = f"output/{base}_y_{fixed_y:.4f}.png"

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_image}")


if __name__ == "__main__":
    main()
