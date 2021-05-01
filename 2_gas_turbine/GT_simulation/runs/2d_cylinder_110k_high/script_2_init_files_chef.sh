#!/bin/bash
# should set up base case to have nominal values of all parameters

ORIG=base_chef
ORIG_TIMESTEP=1
NPROCS=16

for i in {3..200}; do
  D=$(printf 'Run-%03d' $i)

  echo "Setting up run directory $D"
  
  cd "$ORIG" 
  # write geom.smd with correct RVs. 
  ./replaceRVs.sh "$i"

  # Have some pause to ensure file is written over prior to runchef. 
  
  # potentially can remove mesh antics because only .smd is edited. Test this later. 

  # run Chef 
  # only use the geombc.dat....
  cd 16-1-Chef
  ./runChef.sh "$NPROCS"
  # convert to syncio
  ./runConvertO2N.sh
  
  ##/16-procs_case-SyncIO-1/geombc-dat.1

  # Copy generated restart file to respective run directory. 
  cd ../../"$D"/"$NPROCS"-procs_case

  #cp -r ../base_chef/16-1-Chef/16-procs_case .
  cp ../../base_chef/16-1-Chef/16-procs_case-SyncIO-1/geombc-dat.1 .

#  mkdir "$D"/"$NPROCS"-procs_case
#  cd "$D"/"$NPROCS"-procs_case
  

  # return to main directory
  cd ../../
done
