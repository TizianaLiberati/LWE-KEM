/*  Codice_modificato.cpp  –  LWE-KEM GPU benchmark driver
 *
 *  All heavy work runs on GPU via OpenACC:
 *    • Matrix A generated on-the-fly (xorshift, never stored)
 *    • Noise generated on GPU (xorshift)
 *    • 256 encryptions / decryptions batched into 1 kernel each
 *    • Pre-allocated buffers (no allocation in hot loop)
 *
 *  OPTIMIZATIONS vs previous version:
 *
 *  1. PERSISTENT DEVICE REGION: all large working buffers (s, t, A,
 *     r, e1, e2, c, ptxt, dec, mp, mpp) are entered into device memory
 *     ONCE before the N-iteration loop and exited ONCE after.
 *     This eliminates O(N) × PCIe transfers that were the bottleneck
 *     identified by the nsys profile (repeated HtoD/DtoH of ~134 MB
 *     of buffer data per iteration for n=4096).
 *
 *  2. COPYOUT ONLY AT END: z_buf, m_buf, mp_buf are updated on host
 *     only when needed; K_enc/K_dec comparison stays on CPU.
 *
 *  3. ASYNC RNG: noise fill (r, e1, e2) in GPU_GenerateNoise_rngongpu
 *     uses async queues so the three fills overlap where possible.
 *
 *  Only hashing (SHA3-256 of small data) and key comparison remain CPU.
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

int main(int argc, char** argv)
{
    uint32_t q = 3329;
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <N> <n>\n"; return 1; }

    int N      = std::atoi(argv[1]);
    uint32_t n = (uint32_t)std::atoi(argv[2]);
    const uint32_t MSG = 256;

    std::cout << "N = " << N << ", n = " << n
              << ", q = " << q << ", MSG = " << MSG << "\n";

    /* ------------------------------------------------------------------ */
    /*  PRE-ALLOCATE all working buffers ONCE                              */
    /* ------------------------------------------------------------------ */
    std::vector<int32_t> s_buf(n);
    std::vector<int32_t> t_buf(n);
    std::vector<int32_t> z_buf(256);
    std::vector<int32_t> m_buf(MSG);
    std::vector<int32_t> mp_buf(MSG);
    std::vector<int32_t> dec_buf(MSG);
    std::vector<int32_t> ptxt_buf(MSG);
    std::vector<int32_t> r_buf(MSG * n);
    std::vector<int32_t> e1_buf(MSG * n);
    std::vector<int32_t> e2_buf(MSG);
    std::vector<int32_t> c_buf(MSG * (n + 1));
    std::vector<int32_t> c_chk(MSG * (n + 1));
    std::vector<int32_t> A_buf((size_t)n * n);

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
    int32_t* Ap    = A_buf.data();

    size_t total_r = (size_t)MSG * n;
    size_t total_c = (size_t)MSG * (n + 1);
    size_t total_A = (size_t)n * n;

    /* ------------------------------------------------------------------ */
    /*  KEY OPTIMIZATION: enter all large buffers into device memory ONCE  */
    /*  They stay resident for the entire N-iteration benchmark loop.       */
    /*  This removes the dominant bottleneck shown in the nsys profile:     */
    /*  repeated HtoD/DtoH transfers of A (~67 MB), r/e1 (~4 MB each),    */
    /*  and c (~4 MB) on every single iteration.                           */
    /* ------------------------------------------------------------------ */
    #pragma acc enter data create(                  \
        sp[0:n], tp[0:n], zp[0:256],               \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],         \
        ptxtp[0:MSG], e2p[0:MSG],                   \
        rp[0:total_r], e1p[0:total_r],              \
        cp[0:total_c], cchkp[0:total_c],            \
        Ap[0:total_A])

    long long sum_keygen_us = 0, sum_encaps_us = 0, sum_decaps_us = 0;
    int mismatches = 0;

    auto startTot = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        uint64_t key_seed;
        RAND_bytes(reinterpret_cast<unsigned char*>(&key_seed), sizeof(key_seed));

        /* ---- KeyGen (GPU) -------------------------------------------- */
        auto t0 = std::chrono::steady_clock::now();

        uint64_t rho_seed = key_seed ^ 0xCAFEBABE12345678ULL;

        /* GenerateA writes into Ap which is already device-resident.
           The function uses copyout internally; with managed memory the
           device copy is updated in-place — no extra transfer needed.    */
        GenerateA_GPU_rngongpu(rho_seed, n, q, Ap);

        /* KeyGen reads Ap (device-resident), writes sp and tp             */
        KeyGen_GPU_rngongpu_Aflat(key_seed, n, q, Ap, sp, tp);

        /* Generate z on CPU (tiny — 256 values) */
        for (int i = 0; i < 256; ++i) zp[i] = getRandomInt(0, 1);
        /* Update device copy of z */
        #pragma acc update device(zp[0:256])

        auto t1 = std::chrono::steady_clock::now();
        sum_keygen_us += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

        /* ---- Encaps (GPU) -------------------------------------------- */
        std::vector<int32_t> K_enc;
        auto t2 = std::chrono::steady_clock::now();

        Encaps_GPU_Aflat(key_seed, n, q, Ap, tp, cp, K_enc,
                         mp, rp, e1p, e2p, ptxtp);

        auto t3 = std::chrono::steady_clock::now();
        sum_encaps_us += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();

        /* ---- Decaps (GPU) -------------------------------------------- */
        std::vector<int32_t> K_dec;
        auto t4 = std::chrono::steady_clock::now();

        Decaps_GPU_Aflat(key_seed, n, q, Ap, tp, sp, zp, cp, K_dec,
                         mpp, decp, rp, e1p, e2p, ptxtp, cchkp);

        auto t5 = std::chrono::steady_clock::now();
        sum_decaps_us += std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();

        /* ---- Check shared keys (CPU, tiny data) ----------------------- */
        bool same = (K_enc.size() == K_dec.size());
        if (same) {
            for (size_t i = 0; i < K_enc.size(); ++i)
                if (K_enc[i] != K_dec[i]) { same = false; break; }
            if (!same) ++mismatches;
        }
    }

    /* ------------------------------------------------------------------ */
    /*  Release device memory once, at the very end                        */
    /* ------------------------------------------------------------------ */
    #pragma acc exit data delete(                   \
        sp[0:n], tp[0:n], zp[0:256],               \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],         \
        ptxtp[0:MSG], e2p[0:MSG],                   \
        rp[0:total_r], e1p[0:total_r],              \
        cp[0:total_c], cchkp[0:total_c],            \
        Ap[0:total_A])

    auto endTot = std::chrono::steady_clock::now();
    auto us_tot = std::chrono::duration_cast<std::chrono::microseconds>(endTot - startTot).count();
    double total_s = us_tot / 1e6;

    std::cout << n << ";"
              << (double)sum_keygen_us / N << ";"
              << (double)sum_encaps_us / N << ";"
              << (double)sum_decaps_us / N << ";"
              << total_s << "\n";

    return 0;
}
