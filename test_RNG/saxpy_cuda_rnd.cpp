#include "rngongpu/rand_aes/aes_rng.cuh"

#include <cuda_runtime.h>
#include <vector>
#include <cstdint>
#include <cstdio>

#include <openssl/rand.h>

extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma)
{
    //Seleziono AES
    using AesRng = rngongpu::RNG<rngongpu::Mode::AES>;
    static AesRng* drbg = nullptr;

    // Noise e nonce fissato (otteniamo sempre gli stessi valori)
/*     std::vector<std::uint8_t> entropy = {
            0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
            0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,
            0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
            0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F
        };
        
        std::vector<std::uint8_t> nonce = {
            0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27
        }; */
        
        std::vector<std::uint8_t> personalization = {};
 
        if (drbg == nullptr) {
        //Alloco vettore per entropy e nonce
        std::vector<std::uint8_t> entropy(32);
        std::vector<std::uint8_t> nonce(8);
        std::vector<std::uint8_t> personalization = {};
        

        //Sovrascrivo con valori random presi da CPU-OpenSSL (RAND_bytes() genera numeri random crittograficamente sicuri)
        //AES deve avere in input de mueri crittograficamente sicuri e questi li creiamo con RAND_bytes di openssl su CPU
        
        if (RAND_bytes(entropy.data(), (int)entropy.size()) != 1) {
            std::fprintf(stderr, "Errore: RAND_bytes(entropy) fallita\n");
            return;
        }

        if (RAND_bytes(nonce.data(), (int)nonce.size()) != 1) {
            std::fprintf(stderr, "Errore: RAND_bytes(nonce) fallita\n");
            return;
        }
        
        
        drbg = new AesRng(
            entropy,
            nonce,
            personalization,
            rngongpu::SecurityLevel::AES128
        );
    }

    std::vector<std::uint8_t> additional_input = {};

    drbg->normal_random_number(
        sigma,
        reinterpret_cast<f64*>(d_out),
        n,
        additional_input
    );

    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        std::fprintf(stderr, "CUDA error after RNGonGPU call: %s\n",
                     cudaGetErrorString(err));
    }
}