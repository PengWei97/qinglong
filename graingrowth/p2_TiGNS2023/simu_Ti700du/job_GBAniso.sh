#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 2
#SBATCH -n 256

mpiexec -n 256 ~/projects/yinglong/yinglong-opt -i p2_GNS_Ti_AGG_gbAniso_stored_700du.i --recover

#  --recover
# tar -cvf - ex_loc3_v4_1/*.e-s0??[2].* | pigz -9 -p 20 > ex_gbAnisotropyMisori_vt.tgz 
# tar -cvf - csv_noTwin_?0min/* | pigz -9 -p 20 > csv_noTwin.tgz
