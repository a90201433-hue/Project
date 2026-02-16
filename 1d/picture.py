import pandas as pd
import matplotlib.pyplot as plt
import re
import sys 

# --- КОНСТАНТЫ И КОНФИГУРАЦИЯ ---

GAMMA = 1.4
METHOD_STYLES = {
    'Godunov': {'color': 'blue', 'linestyle': '-', 'label': 'Godunov'},
    'Kolgan': {'color': 'red', 'linestyle': '-', 'label': 'Godunov-Kolgan'},
    'ENO': {'color': 'purple', 'linestyle': '-', 'label': 'ENO'},
    'WENO': {'color': 'green', 'linestyle': '-', 'label': 'WENO'},
    'Rodionov': {'color': 'dodgerblue', 'linestyle': '-', 'label': 'Godunov-Kolgan-Rodionov'},
    'Riemann': {'color': 'black', 'linestyle': '-.', 'label': 'Analytical'},
    # Стиль для TVD
    'TVD': {'color': 'navy', 'linestyle': '-', 'label': 'TVD'}
}

# --- ФУНКЦИИ ПАРСИНГА И ЧТЕНИЯ ДАННЫХ ---

def parse_config_scheme(config_path='config.toml'):
    """Читает параметры из секции [scheme] в config.toml, используя более простое регулярное выражение."""
    scheme_params = {
        'methods': [], 
        'TVD': 'Off', 
        'High_order_method': '', 
        'TVD_limiter': ''
    }
    
    try:
        with open(config_path, 'r') as f:
            content = f.read()
            
            # Регулярное выражение для поиска всего блока [scheme]
            scheme_match = re.search(r'\[scheme\]\s*(.*?)(?=\n\[|$)', content, re.DOTALL)
            
            if scheme_match:
                scheme_block = scheme_match.group(1)
                
                # Регулярное выражение для извлечения пар ключ = значение
                line_pattern = re.compile(r'^\s*([a-zA-Z_]+)\s*=\s*(.*?)\s*$', re.MULTILINE)
                
                for key, value_raw in line_pattern.findall(scheme_block):
                    
                    # Очистка значения: удаление внешних кавычек и пробелов
                    value = value_raw.strip().strip('"')
                    
                    if key == 'methods':
                        # Специальный парсинг для массива
                        if value.startswith('[') and value.endswith(']'):
                            value = value.strip('[]')
                            scheme_params['methods'] = [
                                m.strip().strip('"') for m in value.split(',') if m.strip()
                            ]
                        
                    elif key in scheme_params:
                        scheme_params[key] = value
                        
            # 💡 ОТЛАДКА: Выводим, что было прочитано
            
            
    except FileNotFoundError:
        print(f"Ошибка: Файл {config_path} не найден.")
        
    return scheme_params


def read_data_file(name):
    """Считывает данные из CSV-файла и возвращает DataFrame."""
    file_path = "CSV_files/" + name + ".csv"
    try:
        data = pd.read_csv(file_path)
        
        return data
    except FileNotFoundError:
        
        return None
    except pd.errors.EmptyDataError:
        
        return None
    except Exception as e:
       
        return None


def fill_axes(axs, data, method_name, alpha=1.0):
    """Строит графики для заданного набора данных и метода."""
    
    # Доступ к глобальному METHOD_STYLES должен быть здесь
    style_info = METHOD_STYLES.get(method_name, 
                                   {'color': 'gray', 'linestyle': '-', 'label': method_name})
    
    cols = ['rho', 'u', 'P', 'e']
    
    # Проверка столбцов: если столбцы не совпадают, график не строится
    if not all(col in data.columns for col in ['x'] + cols):
       
        return style_info['label'] 
        
    # --- Построение графиков для ВСЕХ ЧЕТЫРЕХ ПЕРЕМЕННЫХ ---
    
    plot_axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
    
    for ax, col in zip(plot_axes, cols):
        ax.plot(data['x'], data[col], 
                color=style_info['color'], 
                linestyle=style_info['linestyle'], 
                alpha=alpha)
        
    # Установка меток для Y-осей
    axs[0, 0].set_ylabel(r'$\rho$')
    axs[0, 1].set_ylabel('u')
    axs[1, 0].set_ylabel('P')
    axs[1, 1].set_ylabel(r'$\varepsilon$')
    
    # Установка меток для X-осей
    for ax in axs.flatten():
        ax.set_xlabel('x')

    return style_info['label']


# --- ОСНОВНАЯ ЛОГИКА ---

def main():
    plt.style.use('seaborn-v0_8')
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True) 

    # 1. Считываем все параметры схемы
    scheme_params = parse_config_scheme()
    methods_list = scheme_params['methods']
    
    # 2. Проверяем флаг TVD и настраиваем метод "TVD"
    if scheme_params.get('TVD', 'Off') == 'On':
        tvd_limiter = scheme_params.get('TVD_limiter', 'DefaultLimiter')
        high_order = scheme_params.get('High_order_method', 'DefaultMethod')
        
        tvd_label = f"TVD with {tvd_limiter} by {high_order}"
        
        
        METHOD_STYLES['TVD']['label'] = tvd_label
        
        # Добавляем "TVD" в список методов
        if 'TVD' not in methods_list:
            methods_list.append('TVD')
    
    # 3. Добавляем аналитическое решение
    analytical_method = 'Riemann'
    if analytical_method not in methods_list:
        methods_list_with_analytical = methods_list + [analytical_method]
    else:
        methods_list_with_analytical = methods_list

    
    
    
    legend_labels = []

    # 4. Цикл построения графиков
    for method_name in methods_list_with_analytical:
        
        alpha = 0.6 if method_name == analytical_method else 1.0
        
        data = read_data_file(method_name)
        
        if data is not None:
            label = fill_axes(axs, data, method_name, alpha)
            legend_labels.append(label)

    # 5. Финальная настройка графиков
    
    fig.legend(legend_labels, loc='center left', bbox_to_anchor=(0.9, 0.5), fontsize=14)

    for ax in axs.flatten():
        ax.grid(True, alpha=0.5)

    plt.tight_layout(rect=[0, 0, 0.9, 1]) 
    plt.savefig(f'output/picture.png', bbox_inches='tight', dpi=300, transparent=False)
    plt.show()

if __name__ == '__main__':
    main()