import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

# === НАСТРОЙКИ ===
BASE_DIR = "./"
METHOD_DIRS = ["CSV_files/WENO"]      # Если оставить [], то методы найдутся автоматически

#METHOD_DIRS = ["CSV_files/Godunov", "CSV_files/Kolgan2", "CSV_files/ENO" ]
X_COLUMN = "x"
VAR_COLUMNS = ["rho", "u", "P", "e"]      # 4 переменные для 2×2 анимации
FPS = 10


# === Автоматический поиск методов ===
if not METHOD_DIRS:
    METHOD_DIRS = [d for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d))]
# print("Методы:", METHOD_DIRS)


# === Считываем данные по методам ===
methods_data = {}

for method in METHOD_DIRS:
    path = os.path.join(BASE_DIR, method)
    files = sorted(
        glob.glob(os.path.join(path, "step_*.csv")),
        key=lambda p: int(re.search(r"step_(\d+)\.csv", p).group(1))
    )

    if not files:
        print(f"В методе {method} нет файлов — пропуск.")
        continue

    dfs = []
    steps = []

    for f in files:
        step = int(re.search(r"step_(\d+)\.csv", f).group(1))
        df = pd.read_csv(f)
        dfs.append(df)
        steps.append(step)

    x = dfs[0][X_COLUMN].to_numpy()

    # создаём словарь переменных → 2D массив (время × координаты)
    var_arrays = {var: np.vstack([df[var].to_numpy() for df in dfs]) for var in VAR_COLUMNS}

    methods_data[method] = {
        "steps": steps,
        "x": x,
        "vars": var_arrays
    }


# Проверка
if not methods_data:
    raise RuntimeError("Нет данных ни для одного метода!")


# === Общий список шагов ===
all_steps = sorted(set().union(*[d["steps"] for d in methods_data.values()]))
nsteps = len(all_steps)
# print("Общие шаги:", all_steps)


# === Подготовка фигуры 2×2 ===
plt.style.use('seaborn-v0_8')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()   # теперь 0..3

titles = [r'$\rho$', "u", "P", "e"]

# линии: dict[var][method] = line2d
lines = {var: {} for var in VAR_COLUMNS}

COLORS = {
    "CSV_files/Godunov": "blue",
    #"CSV_files/Kolgan": "red",
    "CSV_files/Kolgan2": "purple",
    #"CSV_files/Rodionov": "green", 
    "CSV_files/Rodionov2": "cyan",
    "CSV_files/ENO" : "pink",
    "CSV_files/WENO" : "red"

}
LABELS = {
    "CSV_files/Godunov": "Godunov",
    "CSV_files/Kolgan2": "Kolgan2",
    #"CSV_files/Rodionov": "Rodionov",
    "CSV_files/Rodionov2": "Rodionov2",
    "CSV_files/ENO" : "ENO",
    "CSV_files/WENO" : "WENO"
}

# === Инициализация графиков ===
for i, var in enumerate(VAR_COLUMNS):
    ax = axes[i]
    ax.set_title(titles[i])
    # ax.set_xlabel(X_COLUMN)
    # ax.set_ylabel(var)
    ax.grid(True)

    # y-лимиты для переменной — общие по всем методам
    ymin = min(d["vars"][var].min() for d in methods_data.values())
    ymax = max(d["vars"][var].max() for d in methods_data.values())
    dy = ymax - ymin
    ax.set_ylim(ymin - 0.05*dy, ymax + 0.05*dy)

    # создаём линии для всех методов
    for method, data in methods_data.items():
        color = COLORS.get(method, None)
        label = LABELS.get(method, None)
        line, = ax.plot(data["x"], data["vars"][var][0], label=label, c=color)
        lines[var][method] = line

    ax.legend(fontsize=14)


# === Функция обновления ===
def update(frame_i):
    step = all_steps[frame_i]

    for var in VAR_COLUMNS:
        for method, data in methods_data.items():
            if step in data["steps"]:
                idx = data["steps"].index(step)
                lines[var][method].set_ydata(data["vars"][var][idx])
            # если шага нет — линия не обновляется

    fig.suptitle(f"Step {step}", fontsize=16)
    return [line for var in VAR_COLUMNS for line in lines[var].values()]


ani = FuncAnimation(fig, update, frames=nsteps, blit=True, interval=100)


# === Сохранение ===
out_gif = "output/animation.gif"
out_mp4 = "output/animation.mp4"

ani.save(out_gif, writer=PillowWriter(fps=FPS))
print("Сохранено GIF:", out_gif)

try:
    ani.save(out_mp4, writer=FFMpegWriter(fps=FPS))
    print("Сохранено MP4:", out_mp4)
except:
    print("MP4 сохранить не удалось (нет ffmpeg)")


