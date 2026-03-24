#include <iostream>
#include <vector>
#include <openacc.h>
#include <cuda_runtime.h>

#include "rngongpu/rand_aes/aes_rng.cuh"

using f64 = double;

// Forward declaration of CUDA kernel
__global__ void copy_kernel(const f64* src, f64* dst, int N);

int main()
{
    const int N = 65536;

    // --------------------------
    // RNGonGPU setup
    // --------------------------
    std::vector<unsigned char> entropy(32);
    for(int i=0;i<32;i++) entropy[i] = i;

    std::vector<unsigned char> nonce = {0x20,0x21,0x22,0x23,
                                        0x24,0x25,0x26,0x27};

    rngongpu::RNG<rngongpu::Mode::AES> drbg(
        entropy, nonce, {}, rngongpu::SecurityLevel::AES128
    );

    std::cout << "RNG instantiated." << std::endl;
    drbg.print_params();

    // --------------------------
    // Allocate RNGonGPU device array
    // --------------------------
    f64* d_rng = nullptr;
    cudaMalloc(&d_rng, N * sizeof(f64));

    std::vector<unsigned char> additional_input = {};
    drbg.normal_random_number(1.0, d_rng, N, additional_input);

    // --------------------------
    // Allocate OpenACC-managed array
    // --------------------------
    f64* d_acc = nullptr;
    #pragma acc enter data create(d_acc[0:N])

    // --------------------------
    // Copy RNGonGPU array into OpenACC array
    // --------------------------
    f64* d_acc_ptr = (f64*) acc_deviceptr(d_acc);

    int threads = 256;
    int blocks = (N + threads - 1) / threads;

    copy_kernel<<<blocks, threads>>>(d_rng, d_acc_ptr, N);
    cudaDeviceSynchronize();

    // --------------------------
    // OpenACC computation
    // --------------------------
    #pragma acc parallel loop deviceptr(d_acc)
    for(int i=0; i<N; i++)
        d_acc[i] *= 2.0;

    // --------------------------
    // Copy first 10 numbers back to host for demonstration
    // --------------------------
    f64 h[10];
    #pragma acc update self(h[0:10]) async(0)

    std::cout << "Random numbers (first 10, scaled by 2):\n";
    for(int i=0; i<10; i++)
        std::cout << h[i] << "\n";

    // --------------------------
    // Clean up
    // --------------------------
    #pragma acc exit data delete(d_acc)
    cudaFree(d_rng);

    return 0;
}