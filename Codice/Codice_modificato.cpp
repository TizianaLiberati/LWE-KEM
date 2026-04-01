/*  Codice_modificato.cpp  –  LWE-KEM GPU benchmark driver
 *
 *  All heavy work runs on GPU via OpenACC:
 *    • Matrix A generated on-the-fly (xorshift, never stored)
 *    • Noise generated on GPU (xorshift)
 *    • 256 encryptions / decryptions batched into 1 kernel each
 *    • Pre-allocated buffers (no allocation in hot loop)
 *
 *  Only hashing (SHA3-256 of small data) remains on CPU.
 *
 *  Compile:  nvc++ -std=c++20 -O2 -acc -gpu=cc80,managed -Minfo=accel
 *  Run:      ./lwe_kem <N> <n>
 */
#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <random>
#include "pke.h"
#include "kem.h"
#include "utils.h"
#include <openssl/rand.h>

//
//extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma);
//

int main(int argc, char** argv)
{
    uint32_t q = 3329;
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <N> <n>\n"; return 1; }

    int N = std::atoi(argv[1]);
    uint32_t n = (uint32_t)std::atoi(argv[2]);
    const uint32_t MSG = 256;

    //
    const int RNG_TEST_N = 4;
    std::vector<double> rng_test_buf(RNG_TEST_N);
    double* rngp = rng_test_buf.data();
    //

    std::cout << "N = " << N << ", n = " << n
              << ", q = " << q << ", MSG = " << MSG << "\n";

    //
        /* ---- RNGonGPU test ---- */
        /* 
    #pragma acc data copyout(rngp[0:RNG_TEST_N])
    {
        #pragma acc host_data use_device(rngp)
        {
            rngongpu_fill_normal(rngp, RNG_TEST_N, 1.0);
        }
    }

    std::cout << "RNGonGPU test values:\n";
    for (int i = 0; i < RNG_TEST_N; ++i) {
        std::cout << "rng_test_buf[" << i << "] = " << rng_test_buf[i] << "\n";
    } */
    //

    /* ---- PRE-ALLOCATE all working buffers ONCE ---- */
    std::vector<int32_t> s_buf(n);          /* secret key         */
    std::vector<int32_t> t_buf(n);          /* public key         */
    std::vector<int32_t> z_buf(256);        /* implicit rejection  */
    std::vector<int32_t> m_buf(MSG);        /* plaintext bits     */
    std::vector<int32_t> mp_buf(MSG);       /* decrypted bits     */
    std::vector<int32_t> dec_buf(MSG);      /* raw decrypt output */
    std::vector<int32_t> ptxt_buf(MSG);     /* encoded plaintext  */
    std::vector<int32_t> r_buf(MSG * n);    /* noise r            */
    std::vector<int32_t> e1_buf(MSG * n);   /* noise e1           */
    std::vector<int32_t> e2_buf(MSG);       /* noise e2           */
    std::vector<int32_t> c_buf(MSG*(n+1));  /* ciphertext         */
    std::vector<int32_t> c_chk(MSG*(n+1));  /* re-encrypt check   */
    std::vector<int32_t> A_buf((size_t)n * n);

    
    /* raw pointers for GPU (managed memory) */
    int32_t* sp    = s_buf.data();
    int32_t* tp    = t_buf.data();
    int32_t* zp    = z_buf.data();
    int32_t* mp    = m_buf.data();
    int32_t* mpp   = mp_buf.data();
    int32_t* decp  = dec_buf.data();
    int32_t* ptxtp = ptxt_buf.data();
    int32_t* rp    = r_buf.data();
    int32_t* e1p   = e1_buf.data();
    int32_t* e2p   = e2_buf.data();
    int32_t* cp    = c_buf.data();
    int32_t* cchkp = c_chk.data();
    int32_t* Ap = A_buf.data();


    /* CPU RNG for key seeds */
    //std::mt19937_64 seed_rng(42);
    

    long long sum_keygen_us = 0, sum_encaps_us = 0, sum_decaps_us = 0;
    int mismatches = 0;

    auto startTot = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        //uint64_t key_seed = seed_rng();
        uint64_t key_seed;
        RAND_bytes(reinterpret_cast<unsigned char*>(&key_seed), sizeof(key_seed));
        //std::cout << "key_seed = " << key_seed << "\n";
        
        /* ---- KeyGen (GPU) ---- */
        auto t0 = std::chrono::steady_clock::now();
        //KeyGen_GPU(key_seed, n, q, sp, tp);
        //KeyGen_GPU_rngongpu(key_seed, n, q, sp, tp);
        uint64_t rho_seed = key_seed ^ 0xCAFEBABE12345678ULL;
        GenerateA_GPU_rngongpu(rho_seed, n, q, Ap);

        /* if (k == 0) {
            std::cout << "\nFirst values of A_flat:\n";
            for (int idx = 0; idx < 8; ++idx) {
                std::cout << "A_flat[" << idx << "] = " << Ap[idx] << "\n";
            }
        } */
    
    
        KeyGen_GPU_rngongpu_Aflat(key_seed, n, q, Ap, sp, tp);
        /* Generate z on CPU (tiny — 256 values) */
        for (int i = 0; i < 256; ++i) zp[i] = getRandomInt(0, 1);
        auto t1 = std::chrono::steady_clock::now();
        sum_keygen_us += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

        /* ---- Encaps (GPU) ---- */
        std::vector<int32_t> K_enc;
        auto t2 = std::chrono::steady_clock::now();
        //Encaps_GPU(key_seed, n, q, tp, cp, K_enc,
        //           mp, rp, e1p, e2p, ptxtp);
        Encaps_GPU_Aflat(key_seed, n, q, Ap, tp, cp, K_enc, mp, rp, e1p, e2p, ptxtp);
        auto t3 = std::chrono::steady_clock::now();
        sum_encaps_us += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();

        /* ---- Decaps (GPU) ---- */
        std::vector<int32_t> K_dec;
        auto t4 = std::chrono::steady_clock::now();
        //Decaps_GPU(key_seed, n, q, tp, sp, zp, cp, K_dec,
        //          mpp, decp, rp, e1p, e2p, ptxtp, cchkp);
        Decaps_GPU_Aflat(key_seed, n, q, Ap, tp, sp, zp, cp, K_dec,
                   mpp, decp, rp, e1p, e2p, ptxtp, cchkp);
        auto t5 = std::chrono::steady_clock::now();
        sum_decaps_us += std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();

        /* ---- Check shared keys ---- */
        bool same = (K_enc.size() == K_dec.size());
        if (same) {
            for (size_t i = 0; i < K_enc.size(); ++i)
                if (K_enc[i] != K_dec[i]) { same = false; break; }
            if (!same) ++mismatches;
        }
    }

    auto endTot = std::chrono::steady_clock::now();
    auto us_tot = std::chrono::duration_cast<std::chrono::microseconds>(endTot - startTot).count();
    double total_s = us_tot / 1e6;

    std::cout << n << ";"
              << (double)sum_keygen_us / N << ";"
              << (double)sum_encaps_us / N << ";"
              << (double)sum_decaps_us / N << ";"
              << total_s << /* ";"
              << mismatches <<  */"\n";

    return 0;
}
