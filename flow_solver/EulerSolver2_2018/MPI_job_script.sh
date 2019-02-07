#!/bin/bash
###UGEオプション####
#$ -P GR06APR18
#$ -jc dma.LS
#$ -N valid
#$ -cwd
#$ -l h_rt=20:00:00
#$ -pe impi_pslots 560

. /etc/profile.d/modules.sh
module load intel/2018.2.046

mpirun ./bin/Debug/EulerSolver2



