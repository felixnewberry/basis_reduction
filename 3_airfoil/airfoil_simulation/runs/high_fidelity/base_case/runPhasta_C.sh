#/bin/bash

#pconf=/users/skinnerr/git-phasta/buildC-eric-D;   VER=C
#pconf=/users/skinnerr/git-phasta/buildC-master-D; VER=C

pconf=/users/skinnerr/git-phasta/buildC-NACA-R; VER=C
#pconf=/users/skinnerr/git-phasta/buildC-NACA-D; VER=C
#pconf=/users/skinnerr/git-phasta/buildI-NACA-R; VER=IC
#pconf=/users/skinnerr/git-phasta/buildI-NACA-D; VER=IC

rm solver.inp
ln -s solver.inp_C solver.inp

PHASTA_CONFIG=$pconf

date
mpirun -np 16 -x PHASTA_CONFIG=$pconf $pconf/bin/phasta${VER}.exe
date
