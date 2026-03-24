#include <cuda_runtime.h>
using f64 = double;

__global__ void copy_kernel(const f64* src, f64* dst, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < N) dst[idx] = src[idx];
}