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

module load nvhpc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

N=10
n=16384

OUTFILE="out_${SLURM_JOB_ID}_N${N}_n${n}.txt"
ERRFILE="err_${SLURM_JOB_ID}_N${N}_n${n}.txt"

exec > "$OUTFILE" 2> "$ERRFILE"

echo "Threads: $OMP_NUM_THREADS"
echo "Running N=$N n=$n"

date
start_time="$(date -u +%s.%N)"

./lwe_kem "$N" "$n"

exit_code=$?

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"

echo "Exit code: $exit_code"
echo "Total of $elapsed seconds elapsed"

date
