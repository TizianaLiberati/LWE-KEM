#!/bin/bash

# Output file
OUTPUT_FILE="benchmark_results_cpu.csv"
THREADS=72

# Write CSV header
echo "timestamp,N,n,exit_code" > "$OUTPUT_FILE"

export OMP_NUM_THREADS=$THREADS
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Define parameter arrays
# N_values=(10 100)
N_values=(10)
n_values=(32 128 256 512 1024 2048 4096 8192 16384 32768)

# Loop over all combinations
for N in "${N_values[@]}"; do
    for n in "${n_values[@]}"; do
        
        echo "Running ./lwe_kem $N $n (threads=$THREADS)"

        start_time=$(date +%s)

        ./lwe_kem "$N" "$n"
        exit_code=$?

        end_time=$(date +%s)
        runtime=$((end_time - start_time))

        # Save results
        echo "$(date),$N,$n,$exit_code,$runtime" >> "$OUTPUT_FILE"

    done
done

echo "Benchmark completed. Results saved to $OUTPUT_FILE"