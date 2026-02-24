#!/bin/bash

for cfg in config_set/*.toml; do
    name=$(basename "$cfg" .toml)
    ./run_one.sh "$name" || exit 1
done