import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

plt.style.use("seaborn-v0_8-white")


def main():

    if len(sys.argv) < 2:
        print("Использование:")
        print("python3 PlotMap.py <config_name> [field]")
        print("field = rho | P | u | v | speed")
        return

    config_name = sys.argv[1]

    if len(sys.argv) >= 3:
        field_name = sys.argv[2]
    else:
        field_name = "rho"

    base_path = os.path.join("output", config_name)
    csv_path = os.path.join(base_path, "CSV", "Final.csv")
    pics_path = os.path.join(base_path, "pics")

    if not os.path.exists(csv_path):
        print("Файл Final.csv не найден:", csv_path)
        return

    os.makedirs(pics_path, exist_ok=True)

    data = pd.read_csv(csv_path)

    # Берём последнее время
    t_max = data["t"].max()
    data = data[data["t"] == t_max]

    # Уникальные координаты
    x_vals = np.sort(data["x"].unique())
    y_vals = np.sort(data["y"].unique())

    # Функция для разворота поля в 2D
    def reshape_field(name):
        return data.pivot(index="y", columns="x", values=name).values

    rho = reshape_field("rho")
    u   = reshape_field("u")
    v   = reshape_field("v")
    P   = reshape_field("P")
    speed = np.sqrt(u**2 + v**2)

    # Выбор поля
    if field_name == "rho":
        field = rho
        title_field = r"$\rho$"
    elif field_name == "P":
        field = P
        title_field = r"$P$"
    elif field_name == "u":
        field = u
        title_field = r"$u$"
    elif field_name == "v":
        field = v
        title_field = r"$v$"
    elif field_name == "speed":
        field = speed
        title_field = r"$|\vec{u}|$"
    else:
        print("Неизвестное поле.")
        return

    # Для colormap
    X, Y = np.meshgrid(x_vals, y_vals)

    # --- Построение ---
    fig, ax = plt.subplots(figsize=(8, 6))

    cmap = ax.pcolormesh(X, Y, field,
                         shading="auto",
                         cmap="viridis")

    cbar = plt.colorbar(cmap)
    cbar.set_label(title_field)

    # Попробуем streamplot (требует равномерную сетку)
    try:
        ax.streamplot(x_vals, y_vals, u, v,
                      color="white",
                      density=1.2,
                      linewidth=0.8)
    except ValueError:
        # Если сетка неравномерная — fallback на стрелки
        ax.quiver(X[::5, ::5], Y[::5, ::5],
                  u[::5, ::5], v[::5, ::5],
                  color="white",
                  scale=50)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{config_name} | {title_field} | t = {t_max:.4f}")

    plt.tight_layout()

    output_image = os.path.join(pics_path, f"Map_{field_name}.png")
    plt.savefig(output_image, dpi=300)
    plt.close()

    print("Сохранено:", output_image)


if __name__ == "__main__":
    main()
