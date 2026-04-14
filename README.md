# LWE-KEM: Learning With Errors Key Encapsulation Mechanism
Master Thesis project: Tiziana Liberti


A C++ implementation of a Lattice-based Key Encapsulation Mechanism (KEM) using the Learning With Errors (LWE) problem with Fujisaki-Okamoto (FO) transform. This project includes both a baseline implementation and a highly optimized parallelized version.

## Overview

This implementation provides:
- **LWE-based Public Key Encryption (PKE)** - Core lattice encryption scheme
- **Fujisaki-Okamoto Transform** - CCA-secure KEM construction
- **Parallelized Noise Generation** - Efficient discrete Gaussian sampling
- **Optimized Matrix Operations** - Cache-friendly linear algebra


## Project Structure

```
.
├── benchmark.md
├── bench.sh
├── Codice
│   ├── benchmark.sh
│   ├── bench.sh
│   ├── Codice_modificato.cpp
│   ├── command.sbatch
│   ├── hash.cpp
│   ├── hash.h
│   ├── hash_keccak.cpp
│   ├── hash_keccak.h
│   ├── hash_openssl.cpp
│   ├── hash_openssl.h
│   ├── implementation.note
│   ├── KeccakCodePackage
│   ├── kem.cpp
│   ├── kem.h
│   ├── Makefile
│   ├── noise.cpp
│   ├── noise.h
│   ├── pke.cpp
│   ├── pke.h
│   ├── rng.h
│   ├── rng_keccak.h
│   ├── rngongpu_adapter.cu
│   ├── rng_openssl.h
│   ├── sha256.h
│   ├── utils.cpp
│   ├── utils.h
│   └── xorshift.h
├── Codice_modificato.cpp
├── local
├── README.md
├── RNGonGPU
│   ├── benchmark
│   ├── build
│   ├── cmake
│   ├── CMakeLists.txt
│   ├── CONTRIBUTING.md
│   ├── example
│   ├── LICENSE
│   ├── README.md
│   ├── SECURITY.md
│   ├── src
│   ├── test
│   └── thirdparty
├── sha256.h
└── test_RNG
    ├── acc_probe.cpp
    ├── build
    ├── CMakeLists.txt
    ├── Makefile
    ├── openacc_c_main
    ├── openacc_c_main.cpp
    ├── openacc_c_main.o
    ├── saxpy_cuda.cu
    ├── saxpy_cuda.o
    ├── saxpy_cuda_rnd.cpp
    └── TestSimo

```

## Modified Makefile (dual build: with and without NVTX)

- make → no NVTX (default, fastest)
- make nvtx=1 → NVTX enabled
- make profile=off → explicitly disable 


