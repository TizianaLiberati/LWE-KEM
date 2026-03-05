#!/bin/bash

NGPU=1
NCORE=8
NSIMPERGPU=32

# Number of threads per simulation
NTHREAD=$(($NCORE/($NGPU*$NSIMPERGPU)))
if [ $NTHREAD -eq 0 ]; then
    NTHREAD=1
fi
export OMP_NUM_THREADS=$NTHREAD

#############################
# CUDA MPS SETUP
#############################

export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log
mkdir -p $CUDA_MPS_PIPE_DIRECTORY $CUDA_MPS_LOG_DIRECTORY

# Start MPS daemon
nvidia-cuda-mps-control -d

#############################
# Benchmark parameters
#############################

N_values=(10)
n_values=(32 128 256 512 1024)

#############################
# Launch simulations
#############################

for (( i=0; i<$NGPU; i++ )); do

  for (( j=0; j<$NSIMPERGPU; j++ )); do

    (
      id=gpu${i}_sim${j}
      rm -rf $id
      mkdir -p $id
      cd $id

      OUTPUT_FILE="benchmark_results.csv"
      echo "timestamp,N,n,exit_code,runtime" > "$OUTPUT_FILE"

      export CUDA_VISIBLE_DEVICES=$i

      for N in "${N_values[@]}"; do
        for n in "${n_values[@]}"; do

          echo "GPU $i SIM $j : Running ./lwe_kem $N $n"

          start_time=$(date +%s)

          ../lwe_kem "$N" "$n"
          exit_code=$?

          end_time=$(date +%s)
          runtime=$((end_time - start_time))

          echo "$(date),$N,$n,$exit_code,$runtime" >> "$OUTPUT_FILE"

        done
      done

      cd ..
    ) &

  done
done

echo "Waiting for simulations to complete..."
wait

#############################
# Shutdown MPS
#############################

echo quit | nvidia-cuda-mps-control

echo "All benchmarks completed."
