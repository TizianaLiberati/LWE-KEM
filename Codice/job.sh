#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --account=cin_staff
#SBATCH --mem=0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=dcgp_qos_dbg
#SBATCH --time=00:30:00
#SBATCH --job-name=weak_node1

#SBATCH --mail-type=ALL
#SBATCH --mail-user=n.shukla@cineca.it
#SBATCH -o out_%j
#SBATCH -e err_%j

module load nvhpc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "Threads: $OMP_NUM_THREADS"

date
start_time="$(date -u +%s.%N)"

./lwe_kem 10 4096

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed"

date
