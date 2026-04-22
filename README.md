# OpenMP Branch
This branch contains only the OpenMP implementation.
For the GPU version, see the OpenACC branch.

## Build
```
cd Codice
make clean
make
```

## Run
Direct CLI mode:
`./lwe_kem N n`
With:
-	`N` =number of KEM iterations (i.e., number of key generations)
-	`n` = lattice dimension

Configuration mode:
`./lwe_kem config.txt`

The `config.txt` file allows the user to configure parameters. Main configurable parameters include:
- `N` (number of iterations)
- `n` (lattice dimension)
- `q` (modulus)
- noise distributions (`eta`, `sigma`)
CLI mode is intended for quick runs with explicit `N` and `n` values.
The configuration file mode allows the user to define a broader set of parameters, including cryptographic and benchmark settings.

## Example
`./lwe_kem 100 4096`

## Benchmark
`./benchmark.sh`

## Notes
This branch is focused on the OpenMP implementation and CPU-side benchmarking.
For the general project overview, see the `main` branch.
