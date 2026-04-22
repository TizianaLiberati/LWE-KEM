# OpenACC Branch
This branch contains the GPU/OpenACC implementation of the LWE-KEM prototype.

## Build
The RNGonGPU dependency must be built before compiling the main code:
```
cd RNGonGPU
cmake -S . -B build
cmake --build build -j

cd ../Codice
make clean
make
```

## Run
`./lwe_kem N n`

## Benchmark scripts
`./benchmark.sh`  
Runs predefined configurations and stores results in CSV format (`benchmark_results.csv`).

`./bench.sh`  
Produces detailed timing results, including key generation, encapsulation, decapsulation, total execution time, and mismatch count.

## Notes
This branch is focused on the OpenACC implementation and GPU-side benchmarking.
For the general project overview, see the `main` branch.
