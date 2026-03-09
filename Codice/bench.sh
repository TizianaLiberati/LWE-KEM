#!/bin/bash

N_values=(10)
n_values=(32 128 256 512 1024 2048 4096 8192 16384  ) 

RESULTS=()

echo "Running benchmarks..."
echo ""

for N in "${N_values[@]}"; do
  for n in "${n_values[@]}"; do
    echo -n "  N=$N  n=$n ... "
    start=$(date +%s%N)
    ./lwe_kem "$N" "$n"
    exit_code=$?
    end=$(date +%s%N)
    runtime_ms=$(( (end - start) / 1000000 ))
    status=$( [ $exit_code -eq 0 ] && echo "OK" || echo "FAIL" )
    echo "$status (${runtime_ms}ms)"
    RESULTS+=("$N,$n,$runtime_ms,$exit_code,$status")
  done
done

echo ""
echo "┌──────┬───────┬─────────────┬────────┐"
printf "│ %4s │ %5s │ %11s │ %6s │\n" "N" "n" "runtime(ms)" "status"
echo "├──────┼───────┼─────────────┼────────┤"
for row in "${RESULTS[@]}"; do
  IFS=',' read -r N n rt ec status <<< "$row"
  printf "│ %4s │ %5s │ %11s │ %6s │\n" "$N" "$n" "$rt" "$status"
done
echo "└──────┴───────┴─────────────┴────────┘"

total=${#RESULTS[@]}
passed=$(printf '%s\n' "${RESULTS[@]}" | awk -F',' '$5=="OK"' | wc -l)
failed=$(( total - passed ))
echo ""
echo "  Total: $total  |  Passed: $passed  |  Failed: $failed"
