#/bin/bash

#pconf=/users/skinnerr/git-phasta/buildC-eric-D;   VER=C
#pconf=/users/skinnerr/git-phasta/buildC-master-D; VER=C

#pconf=/users/skinnerr/git-phasta/buildC-NACA-R; VER=C
#pconf=/users/skinnerr/git-phasta/buildC-NACA-D; VER=C
pconf=/users/skinnerr/git-phasta/buildI-NACA-R; VER=IC
#pconf=/users/skinnerr/git-phasta/buildI-NACA-D; VER=IC

rm solver.inp
ln -s solver.inp_IC solver.inp

PHASTA_CONFIG=$pconf
NPROCS=16
LOGDIR="./logs"

if [ ! -e "$LOGDIR" ]; then
  echo "Creating new directory for PHASTA logs in $LOGDIR"
  mkdir $LOGDIR 
fi

TS=$(cat $NPROCS-procs_case/numstart.dat)
TS=$((TS)) # Trim whitespace by processing at a number

cp solver.inp $LOGDIR/solver.inp.$TS

PLOG="$LOGDIR/phastalog.$TS"

# overwrite log file first, then append in subsequent calls to tee
date | tee $PLOG
mpirun $1 -np $NPROCS -x PHASTA_CONFIG=$pconf $pconf/bin/phasta${VER}.exe 2>&1 | tee -a $PLOG
date | tee -a $PLOG
