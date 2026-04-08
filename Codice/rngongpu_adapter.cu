#include <cuda_runtime.h>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include "../RNGonGPU/src/include/rngongpu/rand_aes/aes_rng.cuh"

static void build_entropy_and_nonce(uint64_t seed, uint64_t stream_id, std::vector<unsigned char>& entropy, std::vector<unsigned char>& nonce)
{
    //Costruisco entropy e nonce a partire da un seed deterministico
    entropy.resize(32); //32 byte
    nonce.resize(8); //8 byte

    // parto da seed (8 byte). 4 * 8 = 32 byte
    uint64_t x0 = seed;
    uint64_t x1 = seed ^ 0x9E3779B97F4A7C15ULL; //XOR una costante
    uint64_t x2 = seed + 0x517CC1B727220A95ULL; //sommato una costante
    uint64_t x3 = ~seed; //complementare bit a bit

    //costruisco entropy concatenando x0, x1, x2, x3
    std::memcpy(entropy.data() +  0, &x0, 8);
    std::memcpy(entropy.data() +  8, &x1, 8);
    std::memcpy(entropy.data() + 16, &x2, 8);
    std::memcpy(entropy.data() + 24, &x3, 8);

    //riempio con stream_id
    std::memcpy(nonce.data(), &stream_id, 8);
}

extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma, uint64_t seed, uint64_t stream_id)
{
    //Seleziono AES
    using AesRng = rngongpu::RNG<rngongpu::Mode::AES>;

    std::vector<unsigned char> entropy;
    std::vector<unsigned char> nonce;
    build_entropy_and_nonce(seed, stream_id, entropy, nonce);

    std::vector<std::uint8_t> personalization = {}; 
    std::vector<std::uint8_t> additional_input = {};    
    AesRng drbg(
        entropy,
        nonce,
        personalization,
        rngongpu::SecurityLevel::AES128
    );
        
    drbg.normal_random_number(
        sigma,
        reinterpret_cast<f64*>(d_out),
        n,
        additional_input
    );
    cudaError_t err = cudaDeviceSynchronize(); 
    if (err != cudaSuccess) { 
        std::fprintf(stderr, "CUDA error after RNGonGPU call: %s\n", cudaGetErrorString(err)); 
    } 
}

extern "C" void rngongpu_fill_uniform_u32(uint32_t* d_out, int n, uint64_t seed, uint64_t stream_id)
{
    using AesRng = rngongpu::RNG<rngongpu::Mode::AES>;

    std::vector<unsigned char> entropy;
    std::vector<unsigned char> nonce;
    build_entropy_and_nonce(seed, stream_id, entropy, nonce);

    std::vector<std::uint8_t> personalization = {};
    std::vector<std::uint8_t> additional_input = {};

    AesRng drbg(
        entropy,
        nonce,
        personalization,
        rngongpu::SecurityLevel::AES128
    );

    drbg.uniform_random_number(
        d_out,
        n,
        additional_input
    );

    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        std::fprintf(stderr, "CUDA error after RNGonGPU uniform call: %s\n",
                     cudaGetErrorString(err));
    }
}
