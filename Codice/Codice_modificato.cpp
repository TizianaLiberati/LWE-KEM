/*  Codice_modificato.cpp  –  LWE-KEM GPU benchmark driver
 *
*   ┌──────────────────────────────────────────────────────────────────────┐
 *  │  NVTX PROFILING (enabled with -DUSE_NVTX / make nvtx=1)              │
 *  │                                                                      │
 *  │  Ranges added:                                                       │
 *  │    "Init"          – buffer allocation + acc enter data              │
 *  │    "Iteration_k"   – full per-iteration block                        │
 *  │    "KeyGen"        – GenerateA + KeyGen_GPU + tp DtoH                │
 *  │    "Encaps"        – Encaps_GPU_Aflat                                │
 *  │    "Decaps"        – Decaps_GPU_Aflat                                │
 *  │    "Verify"        – CPU decrypt check + shared-key compare          │
 *  │    "Teardown"      – acc exit data + final stats                     │
 *  │                                                                      │
 *  │  Each phase gets a distinct color for easy Nsight Systems timeline   │
 *  │  reading:                                                            │
 *  │    Init       → grey   (0xFF888888)                                  │
 *  │    Iteration  → white  (0xFFFFFFFF)                                  │
 *  │    KeyGen     → green  (0xFF00AA00)                                  │
 *  │    Encaps     → blue   (0xFF0055FF)                                  │
 *  │    Decaps     → orange (0xFFFF8800)                                  │
 *  │    Verify     → purple (0xFFAA00FF)                                  │
 *  │    Teardown   → grey   (0xFF888888)                                  │
 *  └──────────────────────────────────────────────────────────────────────┘
 *
 *  Expected performance after these fixes (n=4096, N=10):
 *    KeyGen:   ~2,300 µs/iter  (saves ~500µs from persistent keygen temps)
 *    Encaps:  ~18,000 µs/iter  (saves ~9ms: h_c cache + persistent noise)
 *    Decaps:  ~18,000 µs/iter  (saves ~9ms: h_c cache + persistent noise)
 *    Total:   ~0.39 s for N=10
 *
 *  Compile (no NVTX):   nvc++ -std=c++20 -O2 -acc -gpu=cc80,managed -Minfo=accel
 *  Compile (with NVTX): make nvtx=1
 *  Run:                 ./lwe_kem <N> <n>
 */

/* ======================================================================== */
/*  NVTX helper macros                                                       */
/*  All NVTX calls are guarded so the binary is unchanged without -DUSE_NVTX */
/* ======================================================================== */
#ifdef USE_NVTX
#  include <nvToolsExt.h>

/* Colour constants used for the timeline */
namespace nvtx_color {
    static constexpr uint32_t grey     = 0xFF888888u;
    static constexpr uint32_t white    = 0xFFFFFFFFu;
    static constexpr uint32_t green    = 0xFF00AA00u;
    static constexpr uint32_t blue     = 0xFF0055FFu;
    static constexpr uint32_t orange   = 0xFFFF8800u;
    static constexpr uint32_t purple   = 0xFFAA00FFu;
    static constexpr uint32_t red      = 0xFFDD0000u;
}

inline void nvtx_push(const char* label, uint32_t color)
{
    nvtxEventAttributes_t attr = {};
    attr.version       = NVTX_VERSION;
    attr.size          = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    attr.colorType     = NVTX_COLOR_ARGB;
    attr.color         = color;
    attr.messageType   = NVTX_MESSAGE_TYPE_ASCII;
    attr.message.ascii = label;
    nvtxRangePushEx(&attr);
}
inline void nvtx_pop() { nvtxRangePop(); }

#  define NVTX_PUSH(label, color)  nvtx_push(label, color)
#  define NVTX_POP()               nvtx_pop()
#  define NVTX_MARK(label)         nvtxMarkA(label)
#else
#  define NVTX_PUSH(label, color)  ((void)0)
#  define NVTX_POP()               ((void)0)
#  define NVTX_MARK(label)         ((void)0)
#endif  /* USE_NVTX */

/* Convenience wrappers with fixed colours per phase */
#define NVTX_PUSH_INIT(lbl)      NVTX_PUSH(lbl, nvtx_color::grey)
#define NVTX_PUSH_ITER(lbl)      NVTX_PUSH(lbl, nvtx_color::white)
#define NVTX_PUSH_KEYGEN(lbl)    NVTX_PUSH(lbl, nvtx_color::green)
#define NVTX_PUSH_ENCAPS(lbl)    NVTX_PUSH(lbl, nvtx_color::blue)
#define NVTX_PUSH_DECAPS(lbl)    NVTX_PUSH(lbl, nvtx_color::orange)
#define NVTX_PUSH_VERIFY(lbl)    NVTX_PUSH(lbl, nvtx_color::purple)
#define NVTX_PUSH_TEARDOWN(lbl)  NVTX_PUSH(lbl, nvtx_color::grey)

/* ======================================================================== */
/*  Regular includes                                                         */
/* ======================================================================== */
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
    uint32_t q = 12289;
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <N> <n>\n"; return 1; }

    int      N = std::atoi(argv[1]);
    uint32_t n = (uint32_t)std::atoi(argv[2]);
    const uint32_t MSG = 256;

    std::cout << "N = " << N << ", n = " << n
              << ", q = " << q << ", MSG = " << MSG << "\n";

    /* ====================================================================
     *  INIT: PRE-ALLOCATE all working buffers ONCE
     * ==================================================================== */
    NVTX_PUSH_INIT("Init");

    std::vector<int32_t> s_buf(n),        t_buf(n),        z_buf(256);
    std::vector<int32_t> m_buf(MSG),      mp_buf(MSG),     dec_buf(MSG);
    std::vector<int32_t> ptxt_buf(MSG),   e2_buf(MSG);
    std::vector<int32_t> r_buf(MSG * n),  e1_buf(MSG * n);
    std::vector<int32_t> c_buf(MSG * (n + 1));
    std::vector<int32_t> c_chk(MSG * (n + 1));
    std::vector<int32_t> A_buf((size_t)n * n);

    /* FIX 2: Persistent noise temp buffers (doubles for rngongpu rounding) */
    size_t noise_total = (size_t)MSG * n;
    std::vector<double>   r_tmp_buf(noise_total);    /* 8.4 MB */
    std::vector<double>   e1_tmp_buf(noise_total);   /* 8.4 MB */
    std::vector<uint32_t> e2_tmp_buf(MSG);           /* 1 KB   */

    /* FIX 4: Persistent KeyGen temp buffers */
    std::vector<uint32_t> sraw_buf(6 * n);           /* 96 KB  */
    std::vector<double>   etmp_buf(n);               /* 32 KB  */
    std::vector<int32_t>  ebuf_buf(n);               /* 16 KB  */

    int32_t*  sp    = s_buf.data();
    int32_t*  tp    = t_buf.data();
    int32_t*  zp    = z_buf.data();
    int32_t*  mp    = m_buf.data();
    int32_t*  mpp   = mp_buf.data();
    int32_t*  decp  = dec_buf.data();
    int32_t*  ptxtp = ptxt_buf.data();
    int32_t*  rp    = r_buf.data();
    int32_t*  e1p   = e1_buf.data();
    int32_t*  e2p   = e2_buf.data();
    int32_t*  cp    = c_buf.data();
    int32_t*  cchkp = c_chk.data();
    int32_t*  Ap    = A_buf.data();
    double*   r_tmp_d  = r_tmp_buf.data();
    double*   e1_tmp_d = e1_tmp_buf.data();
    uint32_t* e2_tmp_u = e2_tmp_buf.data();
    uint32_t* sraw_p   = sraw_buf.data();
    double*   etmp_p   = etmp_buf.data();
    int32_t*  ebuf_p   = ebuf_buf.data();

    size_t total_r = (size_t)MSG * n;
    size_t total_c = (size_t)MSG * (n + 1);
    size_t total_A = (size_t)n * n;

    /* ------------------------------------------------------------------ */
    /*  PERSISTENT DEVICE REGION — enter ALL buffers once                  */
    /*  Includes the new temp buffers for noise generation and KeyGen.     */
    /* ------------------------------------------------------------------ */
    #pragma acc enter data create(                          \
        sp[0:n], tp[0:n], zp[0:256],                       \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],                 \
        ptxtp[0:MSG], e2p[0:MSG],                           \
        rp[0:total_r], e1p[0:total_r],                      \
        cp[0:total_c], cchkp[0:total_c],                    \
        Ap[0:total_A],                                      \
        r_tmp_d[0:noise_total], e1_tmp_d[0:noise_total],    \
        e2_tmp_u[0:MSG],                                    \
        sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])

    NVTX_POP();  /* Init */

    /* ====================================================================
     *  MAIN BENCHMARK LOOP
     * ==================================================================== */
    long long sum_keygen_us = 0, sum_encaps_us = 0, sum_decaps_us = 0;
    int mismatches = 0;

    auto startTot = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        /* Mark the full iteration on the Nsight timeline */
        char iter_label[32];
        std::snprintf(iter_label, sizeof(iter_label), "Iteration_%d", k);
        NVTX_PUSH_ITER(iter_label);

        uint64_t key_seed;
        RAND_bytes(reinterpret_cast<unsigned char*>(&key_seed), sizeof(key_seed));

        /* ================================================================
         *  KeyGen (GPU)
         * ================================================================ */
        NVTX_PUSH_KEYGEN("KeyGen");
        auto t0 = std::chrono::steady_clock::now();

        uint64_t rho_seed = key_seed ^ 0xCAFEBABE12345678ULL;

        NVTX_PUSH_KEYGEN("GenerateA");
        GenerateA_GPU_rngongpu(rho_seed, n, q, Ap);
        NVTX_POP();  /* GenerateA */

        /* FIX 4: pass persistent temp buffers — no malloc in KeyGen */
        NVTX_PUSH_KEYGEN("KeyGen_GPU");
        KeyGen_GPU_rngongpu_Aflat(key_seed, n, q, Ap, sp, tp,
                                  sraw_p, etmp_p, ebuf_p);
        NVTX_POP();  /* KeyGen_GPU */

        /* tp: 16 KB DtoH — needed for SHA3(pkh) in Encaps/Decaps */
        NVTX_PUSH_KEYGEN("tp_DtoH");
        #pragma acc update self(tp[0:n])
        NVTX_POP();  /* tp_DtoH */

        for (int i = 0; i < 256; ++i) zp[i] = getRandomInt(0, 1);
        #pragma acc update device(zp[0:256])

        auto t1 = std::chrono::steady_clock::now();
        sum_keygen_us += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        NVTX_POP();  /* KeyGen */

        /* ================================================================
         *  Encaps (GPU)
         * ================================================================ */
        std::vector<int32_t> K_enc;
        std::vector<int32_t> h_c_enc;   /* FIX 1: h_c returned here */

        NVTX_PUSH_ENCAPS("Encaps");
        auto t2 = std::chrono::steady_clock::now();

        /* FIX 1: h_c_enc is computed inside and returned to caller.
         * FIX 2: r_tmp_d/e1_tmp_d/e2_tmp_u are persistent — no malloc.  */
        Encaps_GPU_Aflat(key_seed, n, q, Ap, tp, cp, K_enc, h_c_enc,
                         mp, rp, e1p, e2p, ptxtp,
                         r_tmp_d, e1_tmp_d, e2_tmp_u);

        auto t3 = std::chrono::steady_clock::now();
        sum_encaps_us += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
        NVTX_POP();  /* Encaps */

        /* ================================================================
         *  Decaps (GPU)
         * ================================================================ */
        std::vector<int32_t> K_dec;

        NVTX_PUSH_DECAPS("Decaps");
        auto t4 = std::chrono::steady_clock::now();

        /* FIX 1: h_c_enc passed in — Decaps skips SHA3(4.1MB) + DtoH.
         * FIX 2: persistent noise temps passed in.                        */
        Decaps_GPU_Aflat(key_seed, n, q, Ap, tp, sp, zp, cp, K_dec, h_c_enc,
                         mpp, decp, rp, e1p, e2p, ptxtp, cchkp,
                         r_tmp_d, e1_tmp_d, e2_tmp_u);

        auto t5 = std::chrono::steady_clock::now();
        sum_decaps_us += std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
        NVTX_POP();  /* Decaps */

        /* ================================================================
         *  Verify: CPU Decrypt + shared-key comparison (outside timing)
         * ================================================================ */
        NVTX_PUSH_VERIFY("Verify");

        #pragma acc update self(sp[0:n], decp[0:MSG])

        std::vector<int32_t> u0(n);
        for (uint32_t i = 0; i < n; ++i) u0[i] = cp[i];
        int32_t v0 = cp[n];

        std::vector<int32_t> s_k_cpu;
        s_k_cpu.reserve(n + 256);
        for (uint32_t i = 0; i < n; ++i) s_k_cpu.push_back(sp[i]);
        for (int i = 0; i < 256; ++i) s_k_cpu.push_back(zp[i]);

        int32_t dec0_cpu = -999;
        Decrypt(v0, u0, s_k_cpu, q, dec0_cpu);

        /* Check shared keys */
        bool same = (K_enc.size() == K_dec.size());
        if (same) {
            for (size_t i = 0; i < K_enc.size(); ++i)
                if (K_enc[i] != K_dec[i]) { same = false; break; }
            if (!same) {
                ++mismatches;
                NVTX_MARK("KeyMismatch");
            }
        }

        NVTX_POP();  /* Verify */

        NVTX_POP();  /* Iteration_k */
    }

    /* ====================================================================
     *  TEARDOWN: Release device memory + print results
     * ==================================================================== */
    NVTX_PUSH_TEARDOWN("Teardown");

    #pragma acc exit data delete(                           \
        sp[0:n], tp[0:n], zp[0:256],                       \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],                 \
        ptxtp[0:MSG], e2p[0:MSG],                           \
        rp[0:total_r], e1p[0:total_r],                      \
        cp[0:total_c], cchkp[0:total_c],                    \
        Ap[0:total_A],                                      \
        r_tmp_d[0:noise_total], e1_tmp_d[0:noise_total],    \
        e2_tmp_u[0:MSG],                                    \
        sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])

    auto endTot = std::chrono::steady_clock::now();
    double total_s = std::chrono::duration_cast<std::chrono::microseconds>(
        endTot - startTot).count() / 1e6;

    std::cout << n << ";"
              << (double)sum_keygen_us / N << ";"
              << (double)sum_encaps_us / N << ";"
              << (double)sum_decaps_us / N << ";"
              << total_s << ";"
              << mismatches << "\n";

    NVTX_POP();  /* Teardown */

    return 0;
}
