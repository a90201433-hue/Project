import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

# === НАСТРОЙКИ ===
BASE_DIR = "./"
CONFIG_PATH = "config.toml"
# *** ИСПРАВЛЕНО: Указан правильный базовый каталог для данных ***
DATA_BASE_DIR = "CSV_files/ActualRes" 
X_COLUMN = "x"
VAR_COLUMNS = ["rho", "u", "P", "e"]
FPS = 10


# === СЛОВАРЬ СТИЛЕЙ (ЦВЕТ + МЕТКА) ===
# Используем цвета Matplotlib 'tab10' для лучшей различимости
METHOD_STYLES = {
    "Godunov": {"color": "blue", "label": "Godunov"},
    "Kolgan": {"color": "red", "label": "Kolgan"},
    "Kolgan2": {"color": "red", "label": "Kolgan 2"},
    "ENO": {"color": "purple", "label": "ENO"},
    "WENO": {"color": "grey", "label": "WENO"},
    "Rodionov": {"color": "green", "label": "Rodionov"},
    "Rodionov2": {"color": "green", "label": "Rodionov 2"},
    # При необходимости добавьте сюда 'Riemann' для аналитического решения, если оно используется в анимации
}


# --- ФУНКЦИЯ ДЛЯ ЧТЕНИЯ КОНФИГА ВРУЧНУЮ ---

def read_config_methods(config_path):
    """Читает список методов из config.toml вручную (без библиотеки toml)."""
    methods_list = []
    try:
        with open(config_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('methods ='):
                    value = line.replace('methods =', '', 1).strip()
                    if value.startswith('[') and value.endswith(']'):
                        value = value.strip('[]')
                        # Разбиваем по запятой, убираем кавычки и пробелы
                        methods_list = [
                            m.strip().strip('"') for m in value.split(',') if m.strip()
                        ]
                    break
    except FileNotFoundError:
        print(f"Ошибка: Файл {config_path} не найден. Использую запасной список.")
        return ["Godunov", "Kolgan2"]
    
    return methods_list if methods_list else ["Godunov", "Kolgan2"]


# --- ГЛАВНАЯ ЛОГИКА ---

# 1. Считываем методы из конфига
METHOD_NAMES = read_config_methods(CONFIG_PATH)
print("Анализируемые методы:", METHOD_NAMES)

# 2. Считываем данные по методам
methods_data = {}

for method_name in METHOD_NAMES:
    # Формируем полный путь к папке метода
    method_dir = os.path.join(DATA_BASE_DIR, method_name)
    path = os.path.join(BASE_DIR, method_dir)
    
    # *** ИСПРАВЛЕНО: Ищем файлы по шаблону 'MethodName_*.csv' ***
    files = sorted(
        glob.glob(os.path.join(path, f"step_*.csv")),
        # *** ИСПРАВЛЕНО: Извлекаем номер шага из конца имени файла (например, _10.csv) ***
        key=lambda p: int(re.search(r"_(\d+)\.csv$", p).group(1))
    )

    if not files:
        print(f"В методе {method_name} (папка {method_dir}) нет файлов — пропуск.")
        continue

    dfs = []
    steps = []

    for f in files:
        # Извлекаем номер шага из конца имени файла
        step = int(re.search(r"_(\d+)\.csv$", f).group(1))
        df = pd.read_csv(f)
        dfs.append(df)
        steps.append(step)

    x = dfs[0][X_COLUMN].to_numpy()

    var_arrays = {var: np.vstack([df[var].to_numpy() for df in dfs]) for var in VAR_COLUMNS}

    methods_data[method_name] = {
        "steps": steps,
        "x": x,
        "vars": var_arrays
    }


# Проверка: если нет данных, выводим ошибку
if not methods_data:
    raise RuntimeError("Нет данных ни для одного метода! Проверьте DATA_BASE_DIR и имена файлов.")


# === Общий список шагов ===
all_steps = sorted(set().union(*[d["steps"] for d in methods_data.values()]))
nsteps = len(all_steps)


# === Подготовка фигуры 2×2 ===
plt.style.use('seaborn-v0_8')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

titles = [r'$\rho$', "u", "P", "e"]
lines = {var: {} for var in VAR_COLUMNS}


# === Инициализация графиков с использованием METHOD_STYLES ===
for i, var in enumerate(VAR_COLUMNS):
    ax = axes[i]
    ax.set_title(titles[i])
    ax.grid(True)
    ax.set_xlabel('x')
    ax.set_ylabel(titles[i]) # Устанавливаем метку для оси Y

    # y-лимиты для переменной — общие по всем методам
    ymin = min(d["vars"][var].min() for d in methods_data.values())
    ymax = max(d["vars"][var].max() for d in methods_data.values())
    dy = ymax - ymin
    ax.set_ylim(ymin - 0.05*dy, ymax + 0.05*dy)
    
    # создаём линии для всех методов
    for method_name, data in methods_data.items():
        # Динамически получаем цвет и метку из METHOD_STYLES
        style = METHOD_STYLES.get(method_name, {"color": "gray", "label": method_name})
        color = style["color"]
        label = style["label"]
        
        # Строим первую точку (t=0)
        line, = ax.plot(data["x"], data["vars"][var][0], label=label, c=color)
        lines[var][method_name] = line

    # *** ИСПРАВЛЕНО: Выносим легенду из осей в общую легенду фигуры ***
    # Легенда будет добавлена в конце, чтобы была только одна общая
    #ax.get_legend().remove()


# === Функция обновления ===
def update(frame_i):
    step = all_steps[frame_i]

    for var in VAR_COLUMNS:
        for method_name, data in methods_data.items():
            if step in data["steps"]:
                idx = data["steps"].index(step)
                lines[var][method_name].set_ydata(data["vars"][var][idx])
            # если шага нет — линия не обновляется

    fig.suptitle(f"Step {step}", fontsize=16)
    return [line for var in VAR_COLUMNS for line in lines[var].values()]


# === Анимация ===
ani = FuncAnimation(fig, update, frames=nsteps, blit=True, interval=100)


# === Легенда (Вынесена вбок) ===
# Получаем все линии и метки с первой оси, чтобы создать общую легенду
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.85, 0.5), fontsize=14)

# Регулировка макета, чтобы освободить место для легенды
plt.tight_layout(rect=[0, 0, 0.8, 1])


# === Сохранение ===
os.makedirs("output", exist_ok=True) 
out_gif = "output/animation.gif"
out_mp4 = "output/animation.mp4"

ani.save(out_gif, writer=PillowWriter(fps=FPS))
print("Сохранено GIF:", out_gif)

try:
    writer_test = FFMpegWriter(fps=FPS)
    ani.save(out_mp4, writer=writer_test)
    print("Сохранено MP4:", out_mp4)
except Exception as e:
    print(f"MP4 сохранить не удалось (возможно, нет ffmpeg). Ошибка: {e}")
