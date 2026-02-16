# Запуск проекта

## 1. Сборка

В корне проекта выполнить:

```bash
make
```

Появится исполняемый файл:

```
./test
```

Если нужно пересобрать с нуля:

```bash
make clean
make
```

---

## 2. Запуск одного конфига

Конфиги лежат в папке `config_set/`.

Запуск:

```bash
./test config_set/<config_name>.toml
```

После расчёта создадутся папки:

```
output/<config_name>/
    ├── CSV/
    │   └── Final.csv
    └── pics/
```

---

## 3. Построение картинок

### 1D срез

```bash
python3 PlotSlice.py <config_name>
```

Срез по конкретному y:

```bash
python3 PlotSlice.py <config_name> <y_value>
```

---

### 2D тепловая карта

```bash
python3 PlotMap.py <config_name> <field>
```

`<field>` = rho | P | u | v | speed

Параметры отрисовки задаются в `plot_config.json`.

Можно задавать:
- `cmap` — цветовая карта
- `mode` — none | stream | quiver
- `density` — плотность линий тока
- `quiver_step` — разрежение стрелок
- `vmin`, `vmax` — пределы цветовой шкалы
- `logscale` — логарифмическая шкала (true/false)

Настройки применяются автоматически при вызове `PlotMap.py`.

---

### Построить все картинки для всех конфигов

Создать файл `make_pics.sh`:

```bash
#!/bin/bash

for cfg in config_set/*.toml; do
    echo "Making pictures for $cfg"

    name=$(basename "$cfg" .toml)

    python3 PlotSlice.py "$name"

    for field in speed P rho u v; do
        python3 PlotMap.py "$name" "$field"
    done
done
```

Сделать исполняемым:

```bash
chmod +x make_pics.sh
```

Запуск:

```bash
./make_pics.sh
```

---

## 4. Запуск расчёта + построение картинок

Создать файл `run_all.sh`:

```bash
#!/bin/bash

for cfg in config_set/*.toml; do
    echo "Running $cfg"

    ./test "$cfg" || exit 1

    name=$(basename "$cfg" .toml)

    python3 PlotSlice.py "$name"

    for field in speed P rho u v; do
        python3 PlotMap.py "$name" "$field"
    done
done
```

Сделать исполняемым:

```bash
chmod +x run_all.sh
```

Запуск:

```bash
./run_all.sh
```

---

## Обычный порядок работы

```bash
make
./run_all.sh
```
