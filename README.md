# LWE-KEM Parallel Implementation (OpenMP & OpenACC)
This repository provides a parallel implementation of a Learning With Errors (LWE)-based Key Encapsulation Mechanism (KEM), targeting both:
-	CPU architectures via OpenMP
-	GPU architectures via OpenACC

The goal is to study performance trade-offs between multicore CPUs and GPU-accelerated designs.
 
## Overview
The project includes:
-	an OpenMP multithreaded implementation for CPUs
-	an OpenACC GPU implementation
-	configurable LWE parameters
-	benchmark scripts for performance evaluation

## Requirements
Tested with:
-	NVIDIA HPC SDK 24.3 (`nvc++`)
-	OpenSSL development libraries (`libssl-dev`)
-	CMake (required only for building RNGonGPU in the OpenACC branch)

For Docker-based GPU execution, the host system must provide:
-	Docker
-	NVIDIA Container Toolkit
-	NVIDIA GPU drivers compatible with the container runtime
 
## Tested on
-	GH200 (Grace Hopper)
-	H100 (AMD EPYC CPU + NVIDIA GPU)
-	A100
-	V100S
 
## Build
### OpenMP version
```bash
cd Codice
make clean
make
```

### OpenACC version
The OpenACC version requires building the RNGonGPU dependency first. If missing, first fetch its submodule (`git submodule update --init --recursive`).

```bash
cd RNGonGPU
cmake -S . -B build
cmake --build build -j

cd ../Codice
make clean
make
```

## Usage
### CLI mode
Run `./lwe_kem N n`
where:
-	`N` = number of iterations
-	`n` = lattice dimension

Example:
`./lwe_kem 100 4096`

### Configuration file mode (available only for OpenMP version):
Run `./lwe_kem config.txt`. The specified configuration file can be used to set parameters such as:
-	`N`
-	`n`
-	modulus `q`
-	noise/error distribution parameters (e.g. `eta`, `sigma`)
 
## Benchmark
### OpenMP
Run `./benchmark.sh`.
This produces a file named `benchmark_results_cpu.csv`
with the following output format:
`N;n;threads;keygen_us;encaps_us;decaps_us;total_s;mismatches`

### OpenACC
The repository includes benchmark helper scripts such as:
```
./benchmark.sh
./bench.sh
```
Depending on the script used, the output may include either:
-	a runtime log over selected (`N`, `n`) values
-	or detailed timing data in the form:
`N;n;keygen_us;encaps_us;decaps_us;total_s;mismatches`
 
## Docker
A multi-stage Dockerfile is provided to compile both the OpenMP and OpenACC implementations.

Build and run the image with:
```bash
docker build -t lwe-kem .
docker run -it --gpus all lwe-kem
```
Omit the `--gpus all` argument for CPU-only OpenMP testing.

Code for both version will be already compiled inside the container under the folders `/workspace/OpenMP` and `/workspace/OpenACC`.
 
## Notes
-	The OpenACC implementation uses GPU offloading and builds the RNGonGPU dependency separately.
-	The OpenMP and OpenACC versions are compiled with the same toolchain (`nvc++`) for fair comparison.
-	Performance and correctness depend on the selected parameter set, including `N`, `n`, `q`, and the chosen noise parameters/distributions.
