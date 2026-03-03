#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=36
#SBATCH -A eupuser
#SBATCH --nodelist=ngnode02
#SBATCH --partition=t-grace-hopper
#SBATCH --time=24:00:00
module load gcc/gcc-13.3.0
module load nvhpc 
make clean
make
make bench