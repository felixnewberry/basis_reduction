#!/bin/bash
# should set up base case to have nominal values of all parameters

ORIG=base_case
ORIG_TIMESTEP=500
NPROCS=16

for i in {4..200}; do
  D=$(printf 'Run-%03d' $i)

  echo "Setting up run directory $D"

  cp "$ORIG"/solver.inp "$D"

done
