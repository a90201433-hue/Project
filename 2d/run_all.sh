#!/bin/bash

for cfg in config_set/*.toml; do
    echo "Running $cfg"
    echo \n

    ./test "$cfg" || exit 1

    name=$(basename "$cfg" .toml)
    python3 PlotSlice.py "$name"
done
