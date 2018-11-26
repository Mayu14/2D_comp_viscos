#!/bin/bash
###UGEオプション####
#$ -P GR06APR18
#$ -jc dma.A
#$ -N pretest
#$ -cwd
#$ -l h_rt=99:00:00
#$ -pe impi_pslots 40

. /etc/profile.d/modules.sh
module load intel/2018.2.046

mpirun ./bin/Debug/EulerSolver2



