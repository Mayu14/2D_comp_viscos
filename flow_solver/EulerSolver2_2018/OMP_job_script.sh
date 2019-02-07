#!/bin/bash
###UGEオプション####
#$ -P GR06APR18
#$ -jc smb.A
#$ -N pretest
#$ -cwd
#$ -l h_rt=99:00:00
#$ -pe OpenMP 80

. /etc/profile.d/modules.sh
module load intel/2018.2.046

./bin/Debug/EulerSolver2



