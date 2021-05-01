#!/bin/bash
# should set up base case to have nominal values of all parameters

ORIG=base_case
NPROCS=8

for i in {1..5}; do
  D=$(printf 'Run-%03d' $i)

  echo "Setting up run directory $D"

  cp "$ORIG"/solver.inp "$D"

done
