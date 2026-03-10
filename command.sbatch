#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=36
#SBATCH --gres=gpu:1
#SBATCH -A eupuser
#SBATCH --nodelist=ngnode02
#SBATCH --partition=t-grace-hopper
#SBATCH --time=24:00:00

module load gcc/gcc-13.3.0
module load nvhpc

echo "=== Building OpenMP Target Offloading version ==="
make clean
make

echo ""
echo "=== Benchmark (OpenMP target) ==="
make bench
