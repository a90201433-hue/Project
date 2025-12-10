import pandas as pd
import matplotlib.pyplot as plt

# Встроенная библиотека для работы с конфигурационными файлами (INI-подобный формат)
# Используется ручной парсинг, так как нет библиотеки toml

# --- КОНСТАНТЫ И КОНФИГУРАЦИЯ ---

GAMMA = 1.4

# СЛОВАРЬ ДЛЯ АВТОМАТИЧЕСКОГО НАЗНАЧЕНИЯ ЦВЕТОВ И СТИЛЕЙ
# Если нужно больше методов, просто добавьте их сюда.
METHOD_STYLES = {
    'Godunov': {'color': 'blue', 'linestyle': '-', 'label': 'Godunov'},
    'Kolgan': {'color': 'red', 'linestyle': '-', 'label': 'Godunov-Kolgan'},
    'ENO': {'color': 'purple', 'linestyle': '-', 'label': 'ENO'},
    'WENO': {'color': 'green', 'linestyle': '-', 'label': 'WENO'},
    'Rodionov': {'color': 'cyan', 'linestyle': '--', 'label': 'Godunov-Kolgan-Rodionov'},
    # Аналитическое решение (для него устанавливаем alpha=0.6)
    'Riemann': {'color': 'black', 'linestyle': '-.', 'label': 'Analytical'}
}

# --- ГЛОБАЛЬНЫЕ ДАННЫЕ (ОЧИЩЕНО) ---

# Вместо 'global' в функциях будем возвращать данные
def read_config_methods(config_path='config.toml'):
    """Читает список методов из config.toml вручную."""
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
        print(f"Ошибка: Файл {config_path} не найден.")
        # Запасной список
        return ['TVD']
    
    # Если список в конфиге пуст, используем запасной
    return methods_list if methods_list else ['TVD']


def read_data_file(name):
    """Считывает данные из CSV-файла и возвращает DataFrame."""
    file_path = "CSV_files/" + name + ".csv"
    try:
        data = pd.read_csv(file_path)
        return data
    except FileNotFoundError:
        print(f"Предупреждение: Не найден файл CSV для метода '{name}' по пути '{file_path}'.")
        return None


def fill_axes(axs, data, method_name, alpha=1.0):
    """Строит графики для заданного набора данных и метода."""
    
    style_info = METHOD_STYLES.get(method_name, 
                                   {'color': 'gray', 'linestyle': '-', 'label': method_name})
    
    # --- ВАЖНОЕ ИЗМЕНЕНИЕ: используем именованные аргументы color и linestyle ---
    
    # Плотность (rho)
    axs[0, 0].plot(data['x'], data['rho'], 
                   color=style_info['color'], 
                   linestyle=style_info['linestyle'], 
                   alpha=alpha)
    axs[0, 0].set_xlabel('x')
    axs[0, 0].set_ylabel(r'$\rho$')
    
    # Скорость (u)
    axs[0, 1].plot(data['x'], data['u'], 
                   color=style_info['color'], 
                   linestyle=style_info['linestyle'], 
                   alpha=alpha)
    axs[0, 1].set_xlabel('x')
    axs[0, 1].set_ylabel('u')

    # Давление (P)
    axs[1, 0].plot(data['x'], data['P'], 
                   color=style_info['color'], 
                   linestyle=style_info['linestyle'], 
                   alpha=alpha)
    axs[1, 0].set_xlabel('x')
    axs[1, 0].set_ylabel('P')

    # Энергия (e)
    axs[1, 1].plot(data['x'], data['e'], 
                   color=style_info['color'], 
                   linestyle=style_info['linestyle'], 
                   alpha=alpha)
    axs[1, 1].set_xlabel('x')
    axs[1, 1].set_ylabel(r'$\varepsilon$')
    
    return style_info['label']


# --- ОСНОВНАЯ ЛОГИКА ---

def main():
    plt.style.use('seaborn-v0_8')
    fig, axs = plt.subplots(2, 2, figsize=(9, 7), sharex=True)

    # 1. Считываем методы из конфига
    methods_list = read_config_methods()

    # 2. Добавляем аналитическое решение в конец списка
    analytical_method = 'Riemann'
    if analytical_method not in methods_list:
        methods_list_with_analytical = methods_list + [analytical_method]
    else:
        methods_list_with_analytical = methods_list

    legend_labels = []

    # 3. Цикл построения графиков
    for method_name in methods_list_with_analytical:
        
        # Определяем прозрачность (alpha)
        alpha = 0.6 if method_name == analytical_method else 1.0
        
        data = read_data_file(method_name)
        
        if data is not None:
            label = fill_axes(axs, data, method_name, alpha)
            legend_labels.append(label)

    # 4. Финальная настройка графиков
    
    # Единственная легенда вне графиков
    fig.legend(legend_labels, loc='center left', bbox_to_anchor=(0.9, 0.5), fontsize=12)

    for ax in axs.flatten():
        ax.grid(True, alpha=0.5)

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f'output/picture.png', bbox_inches='tight', dpi=300, transparent=False)
    plt.show()

if __name__ == '__main__':
    main()
