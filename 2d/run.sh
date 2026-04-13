#!/bin/bash

EXEC=./test
CONFIG=config_set/config_G.toml

for p in {1..10}
do
    echo "Running with $p processes"
    mpirun -np $p ./test $CONFIG

done