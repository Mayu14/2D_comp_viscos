#!/bin/bash
###UGEオプション####
#$ -P GR06APR18
#$ -jc smb.A
#$ -N 20_800_0010_0200
#$ -cwd
#$ -l h_rt=99:00:00
#$ -pe OpenMP 40

. /etc/profile.d/modules.sh
module load intel/2018.2.046

./ES2_20_800_0010_0200



