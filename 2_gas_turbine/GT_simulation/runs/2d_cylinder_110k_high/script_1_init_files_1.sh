#!/bin/bash
# should set up base case to have nominal values of all parameters

ORIG=base_case
ORIG_TIMESTEP=500
NPROCS=16

for i in {1..200}; do
  D=$(printf 'Run-%03d' $i)

  echo "Setting up run directory $D"

  cp "$ORIG"/solver.inp "$D"
  ##cp "$ORIG"/flowPosix.pht "$D"
  cp "$ORIG"/flowPosix_extract_1.phts "$D"
  cp "$ORIG"/runPhasta.sh "$D"

  
  mkdir "$D"/"$NPROCS"-procs_case
  cd "$D"/"$NPROCS"-procs_case
  
  # Start of a new block; fork from base case
  cp    ../../"$ORIG"/"$NPROCS"-procs_case/numstart.dat .
  cp    ../../"$ORIG"/"$NPROCS"-procs_case/numpe.in .
  #ln -s ../../"$ORIG"/"$NPROCS"-procs_case/geombc.dat.* .
  ln -s ../../"$ORIG"/"$NPROCS"-procs_case/restart-dat."$ORIG_TIMESTEP".* .
  
  #rm -r "$D"/"$NPROCS"-procs_case  
  cd ../../

done
