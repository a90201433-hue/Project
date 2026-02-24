import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import json

plt.style.use("seaborn-v0_8-white")


def main():

    if len(sys.argv) < 3:
        print("Использование:")
        print("python3 PlotMap.py <config_name> <field> [csv_name]")
        return

    config_name = sys.argv[1]
    field_name = sys.argv[2]

    # --- CSV по умолчанию ---
    csv_name = "Final.csv"
    if len(sys.argv) >= 4:
        csv_name = sys.argv[3]

    # --- Читаем конфиг визуализации ---
    with open("plot_config.json", "r") as f:
        plot_cfg = json.load(f)

    if field_name not in plot_cfg:
        print(f"Поле '{field_name}' не найдено в plot_config.json")
        return

    field_cfg = plot_cfg[field_name]

    cmap_name = field_cfg.get("cmap", "viridis")
    mode = field_cfg.get("mode", "none")
    density = field_cfg.get("density", 1.0)
    quiver_step = field_cfg.get("quiver_step", 5)

    base_path = os.path.join("output", config_name)
    csv_path = os.path.join(base_path, "CSV", csv_name)
    pics_path = os.path.join(base_path, "pics")

    if not os.path.exists(csv_path):
        print("CSV файл не найден:", csv_path)
        return

    os.makedirs(pics_path, exist_ok=True)

    data = pd.read_csv(csv_path)

    t_max = data["t"].max()
    data = data[data["t"] == t_max]

    x_vals = np.sort(data["x"].unique())
    y_vals = np.sort(data["y"].unique())

    def reshape_field(name):
        return data.pivot(index="y", columns="x", values=name).values

    rho = reshape_field("rho")
    u   = reshape_field("u")
    v   = reshape_field("v")
    P   = reshape_field("P")
    speed = np.sqrt(u**2 + v**2)

    fields = {
        "rho": (rho, r"$\rho$"),
        "P": (P, r"$P$"),
        "u": (u, r"$u$"),
        "v": (v, r"$v$"),
        "speed": (speed, r"$|\vec{u}|$")
    }

    if field_name not in fields:
        print("Неизвестное поле.")
        return

    field, title_field = fields[field_name]

    X, Y = np.meshgrid(x_vals, y_vals)

    fig, ax = plt.subplots(figsize=(8, 6))

    cmap = ax.pcolormesh(X, Y, field,
                         shading="auto",
                         cmap=cmap_name)

    cbar = plt.colorbar(cmap)
    cbar.set_label(title_field)

    # --- Отображение скорости ---
    if mode == "stream":
        try:
            ax.streamplot(x_vals, y_vals, u, v,
                          color="white",
                          density=density,
                          linewidth=0.8)
        except ValueError:
            ax.quiver(X[::quiver_step, ::quiver_step],
                      Y[::quiver_step, ::quiver_step],
                      u[::quiver_step, ::quiver_step],
                      v[::quiver_step, ::quiver_step],
                      color="white",
                      scale=50)

    elif mode == "quiver":
        ax.quiver(X[::quiver_step, ::quiver_step],
                  Y[::quiver_step, ::quiver_step],
                  u[::quiver_step, ::quiver_step],
                  v[::quiver_step, ::quiver_step],
                  color="white",
                  scale=50)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{config_name} | {title_field} | t = {t_max:.4f}")

    plt.tight_layout()

    # имя картинки = имя csv + поле
    image_name = f"Map_{os.path.splitext(csv_name)[0]}_{field_name}.png"
    output_image = os.path.join(pics_path, image_name)

    plt.savefig(output_image, dpi=300)
    plt.close()

    print("Сохранено:", output_image)


if __name__ == "__main__":
    main()