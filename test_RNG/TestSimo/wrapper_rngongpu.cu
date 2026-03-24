//#include "rngongpu/rand_aes/aes_rng.cuh"
#include <vector>
#include <cstdint>

extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma)
{
    /*
    using AesRng = rngongpu::RNG<rngongpu::Mode::AES>;

    static AesRng* drbg = nullptr;

    if (drbg == nullptr) {
        std::vector<std::uint8_t> entropy = {
            0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
            0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,
            0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
            0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F
        };

        std::vector<std::uint8_t> nonce = {
            0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,
            0x28,0x29,0x2A,0x2B,0x2C,0x2D,0x2E,0x2F
        };

        std::vector<std::uint8_t> personalization = {};

        drbg = new AesRng(
            entropy,
            nonce,
            personalization,
            rngongpu::SecurityLevel::AES128
        );
    }
    std::vector<std::uint8_t> additional_input = {};
    drbg->normal_random_number(sigma, d_out, n, additional_input);
    */
    
    std::cout << "CIAOO2" << std::endl;

    //TODO
    //CREA ARRAY SU GPU E RIEMPIO CON ZERI E POI LO CHIAMO DA ACC_PROBE1.CPP E INCLUDE LA LIBRERIA RNGONGPU
}
