#!/bin/bash

# Проверка аргумента
if [ -z "$1" ]; then
    echo "Использование:"
    echo "./run_one.sh config_name"
    echo "или"
    echo "./run_one.sh config_name.toml"
    exit 1
fi

# Если передали имя с .toml — убираем расширение
cfg_input="$1"
name=$(basename "$cfg_input" .toml)

cfg_path="config_set/$name.toml"

# Проверяем существование файла
if [ ! -f "$cfg_path" ]; then
    echo "Файл $cfg_path не найден"
    exit 1
fi
echo "Running for $name"
echo ""

# Запуск расчёта
./test "$cfg_path" || exit 1

echo ""
echo "Making slice plot"
python3 PlotSlice.py "$name"

echo "Making maps"
for field in speed P rho u v; do
    python3 PlotMap.py "$name" "$field"
done
echo ""
