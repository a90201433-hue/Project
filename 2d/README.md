# Как запускать проект

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

После расчёта создаётся папка:

```
output/<config_name>/
    ├── CSV/
    │   └── Final.csv
    └── pics/
```

---

## 3. Построение графика (если расчёт уже сделан)

```bash
python3 PlotSlice.py <config_name>
```

График сохраняется в:

```
output/<config_name>/pics/Final.png
```

Если нужен срез по конкретному y:

```bash
python3 PlotSlice.py <config_name> 0.5
```

---

## 4. Запуск всех конфигов

Создать файл `run_all.sh`:

```bash
#!/bin/bash

for cfg in config_set/*.toml; do
    ./test "$cfg" || exit 1
    name=$(basename "$cfg" .toml)
    python3 PlotSlice.py "$name"
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

## Обычный порядок работы

```bash
make
./run_all.sh
```
