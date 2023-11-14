#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 2
#SBATCH -n 256

mpiexec -n 256 ~/projects/yinglong/yinglong-opt -i p2_GNS_Ti_AGG_gbAniso_stored_700du.i --recover
