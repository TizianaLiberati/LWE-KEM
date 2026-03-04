#!/bin/bash
for i in $(seq 1 10);
do
    #start_time=$(date +%s%3N)  # start time in milliseconds
    start_time=$(date +%s%6N)  # start time in microseconds

    ./kem  # <-- code for measuring is here

    #end_time=$(date +%s%3N)  # # end time in milliseconds
    end_time=$(date +%s%6N)  # # end time in microseconds

    #duration_ms=$((end_time - start_time))  # duration in milliseconds
    duration_mus=$((end_time - start_time))  # duration in microseconds

    #echo "Execution time in ms: $duration_ms"
    echo "Execution time in mus: $duration_mus"

done