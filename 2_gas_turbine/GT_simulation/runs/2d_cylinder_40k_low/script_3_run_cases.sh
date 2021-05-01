#!/bin/bash

function preComment {
  echo "*****************************************************************************"
  echo "*****************************************************************************"
  echo "*** $1"
}

function echoComment {
  echo "*****************************************************************************"
  echo "*****************************************************************************"
  echo "*** $1"
  echo "*** $(date)"
  echo "*****************************************************************************"
  echo "*****************************************************************************"
}

NPROCS=8

#if (("$NPROCS" < 2)); then
#  MPI_OPTIONS="--bind-to none --mca btl self -np 1"
#else
#  MPI_OPTIONS="-np $NPROCS"
#fi

# Sanitize input
JUST_RUN_CURRENT_CASES_LONGER=true
DIR_INCREMENT=1
NUM_START=1
NUM_RUN_DIRS=200
if [ -z "$1" ] || (("$1" < 1)) || (("$1" > "$DIR_INCREMENT")); then
  echo "Need an input between 1 and $DIR_INCREMENT; this is the starting offset for which directory to run."
  exit
else
  DIR_OFFSET="$1"
fi

if [ "$JUST_RUN_CURRENT_CASES_LONGER" = false ]; then
  echo "Are you sure you want to fork runs from previous runs?"
  read -n1 -r -p "Press y to continue, or any other key to abort..." key
  if [ ! "$key" = 'y' ]; then
    echo
    echo 'Exiting'
    exit
  fi
fi
echo 'Continuing with forking runs from previous runs!'

STARTDIR=$(pwd)

RAND_PARAM_FILE="setting_param_file.dat"
SOLVERINP="solver.inp"

#
# Loop over run directories
#
for (( i=$(($NUM_START+$DIR_OFFSET-1)); i<$(($NUM_START+$NUM_RUN_DIRS)); i+="$DIR_INCREMENT" )); do

  #
  # Read certain settings from file
  #
  #S_TIME_STEPS=$(cat "setting_time_steps.dat")
  #S_TS_SIZE=$(cat "setting_time_step_size.dat")

  ## Replace with matlab entires
  #S_Time_STEPS="$t_step" 
  #S_TS_SIZE="$num_steps"
  # Dive into run directory
  # (***FORMAT SPECIFIES LEADING ZEROS IN Run-0001 etc. DIRECTORIES***)
  #
  RUNDIR=$(printf "Run-%03i" $i)
  cd "$RUNDIR"

  echoComment "Running IC solver for 2D-APGD in $RUNDIR (directory offset is $DIR_OFFSET)"


  ./runPhasta.sh

  #
  # Pop up to starting directory
  #
  cd "$STARTDIR"

done
