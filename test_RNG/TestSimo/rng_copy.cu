#include <cuda_runtime.h>
using f64 = double;

__global__ void copy_kernel(const f64* src, f64* dst, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < N) dst[idx] = src[idx];
}

// Helper function called from main.cpp
void copy_rng_to_acc(const f64* d_rng, f64* d_acc, int N) {
    // Get raw device pointer of OpenACC array
    f64* d_acc_ptr = (f64*) acc_deviceptr(d_acc);

    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    copy_kernel<<<blocks, threads>>>(d_rng, d_acc_ptr, N);
    cudaDeviceSynchronize();
}