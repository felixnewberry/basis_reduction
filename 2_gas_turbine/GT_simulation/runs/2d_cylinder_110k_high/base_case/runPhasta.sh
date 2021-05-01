#/bin/bash

#pconf=/users/fnewberry/phasta/buildI-master-next-ICduct-D; VER=IC
pconf=/users/fnewberry/phasta/buildI-master-next-ICduct-R; VER=IC

#pconf=/users/fnewberry/phasta/buildI-master-next-R; VER=IC
#pconf=/users/skinnerr/git-phasta/buildI-master-next-R; VER=IC
#pconf=/users/skinnerr/git-phasta/buildI-master-next-D; VER=IC

#pconf=/users/skinnerr/git-phasta/buildC-master-next-D; VER=C
#pconf=/users/skinnerr/git-phasta/buildC-master-next-R; VER=C
#pconf=/users/skinnerr/git-phasta/buildC-master-D; VER=C
#pconf=/users/skinnerr/git-phasta/buildC-master-R; VER=C
#pconf=/users/skinnerr/git-phasta/buildC-master-24b5caa0b9bee-R/; VER=C

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
mpirun $1 --bind-to none -np $NPROCS -x PHASTA_CONFIG=$pconf $pconf/bin/phasta${VER}.exe 2>&1 | tee -a $PLOG
date | tee -a $PLOG
