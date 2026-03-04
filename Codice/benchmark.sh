#!/bin/bash

# Output file
OUTPUT_FILE="benchmark_results.csv"

# Write CSV header
echo "timestamp,N,n,exit_code" > "$OUTPUT_FILE"

# Define parameter arrays
N_values=(10 100)
n_values=(32 128 256 512 1024 2048)

# Loop over all combinations
for N in "${N_values[@]}"; do
    for n in "${n_values[@]}"; do
        
        echo "Running ./lwe_kem $N $n"

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
