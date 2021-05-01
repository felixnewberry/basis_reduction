#!/bin/bash

#mpirun -bynode --prefix /usr/local/openmpi/1.6.5-gnu482-thread --hostfile ~/hostfile16 -x PATH -x LD_LIBRARY_PATH -x PATH -x SIM_LICENSE_FILE -np 8 /users/mrasquin/SCOREC.develop-2014/CMake.phParAdapt/build_GNU_OptG/phParAdapt  2>&1 | tee phParAdapt.log

mpirun -np $1 /projects/tools/SCOREC-core/build-viz003/test/chef 2>&1 | tee chef.log


