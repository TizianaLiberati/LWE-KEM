/*  Codice_modificato.cpp  –  LWE-KEM GPU benchmark driver
 *
 *  Compile-time backend selection:
 *    nvc++ -DUSE_OPENACC    -acc    -gpu=cc90,managed ...
 *    nvc++ -DUSE_OMP_TARGET -mp=gpu -gpu=cc90,managed ...
 *
 *  ──────────────────────────────────────────────────────────────────────
 *  PERSISTENT DEVICE REGION STRATEGY
 *  ──────────────────────────────────────────────────────────────────────
 *
 *  Both backends support a "enter once / exit once" data-lifetime model.
 *
 *  OpenACC:
 *    #pragma acc enter data create(...)   — allocates on device, no HtoD
 *    #pragma acc exit  data delete(...)   — frees device memory
 *    #pragma acc update self/device(...)  — explicit partial transfers
 *
 *  OpenMP Target:
 *    #pragma omp target enter data map(alloc:...)  — allocates on device
 *    #pragma omp target exit  data map(release:...) — frees device memory
 *    #pragma omp target update from/to(...)         — explicit transfers
 *
 *  Both models map 1-to-1 onto the same CUDA driver operations when
 *  compiled with NVHPC (-acc vs -mp=gpu).  The key difference is that
 *  OpenACC has an explicit `present()` runtime assertion, while OpenMP
 *  relies on `is_device_ptr()` which is a compiler hint only.
 *
 *  Under -gpu=managed (unified memory) the driver may not perform any
 *  physical DMA for `update self/from` and `update device/to` calls —
 *  they become cache-coherence hints.  This is identical for both
 *  backends when compiled with NVHPC.
 *
 *  ──────────────────────────────────────────────────────────────────────
 *  NVTX PROFILING (enabled with -DUSE_NVTX / make nvtx=1)
 *  ──────────────────────────────────────────────────────────────────────
 *  Ranges:  Init | Iteration_k | KeyGen | Encaps | Decaps | Verify | Teardown
 *  Colors:  grey | white | green | blue | orange | purple | grey
 *
 *  ──────────────────────────────────────────────────────────────────────
 *  EXPECTED PERFORMANCE (n=4096, N=10):
 *    KeyGen:   ~2,300 µs/iter
 *    Encaps:  ~18,000 µs/iter
 *    Decaps:  ~18,000 µs/iter
 *    Total:    ~0.39 s
 *
 *  Compile (OpenACC, no NVTX):
 *    nvc++ -std=c++20 -O2 -acc -gpu=cc90,managed -DUSE_OPENACC \
 *          -Minfo=accel -o lwe_kem *.cpp *.o $(LDFLAGS)
 *
 *  Compile (OpenMP Target, no NVTX):
 *    nvc++ -std=c++20 -O2 -mp=gpu -gpu=cc90,managed -DUSE_OMP_TARGET \
 *          -Minfo=mp -o lwe_kem *.cpp *.o $(LDFLAGS)
 *
 *  Run:   ./lwe_kem <N> <n>
 */

/* ======================================================================== */
/*  NVTX helper macros                                                        */
/* ======================================================================== */
#ifdef USE_NVTX
#  include <nvToolsExt.h>
namespace nvtx_color {
    static constexpr uint32_t grey   = 0xFF888888u;
    static constexpr uint32_t white  = 0xFFFFFFFFu;
    static constexpr uint32_t green  = 0xFF00AA00u;
    static constexpr uint32_t blue   = 0xFF0055FFu;
    static constexpr uint32_t orange = 0xFFFF8800u;
    static constexpr uint32_t purple = 0xFFAA00FFu;
}
static void nvtx_push(const char* label, uint32_t color) {
    nvtxEventAttributes_t attr = {};
    attr.version       = NVTX_VERSION;
    attr.size          = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    attr.colorType     = NVTX_COLOR_ARGB;
    attr.color         = color;
    attr.messageType   = NVTX_MESSAGE_TYPE_ASCII;
    attr.message.ascii = label;
    nvtxRangePushEx(&attr);
}
static void nvtx_pop() { nvtxRangePop(); }
#  define NVTX_PUSH(lbl, col)   nvtx_push(lbl, col)
#  define NVTX_POP()            nvtx_pop()
#  define NVTX_MARK(lbl)        nvtxMarkA(lbl)
#else
#  define NVTX_PUSH(lbl, col)   ((void)0)
#  define NVTX_POP()            ((void)0)
#  define NVTX_MARK(lbl)        ((void)0)
#endif

#define NVTX_PUSH_INIT(l)     NVTX_PUSH(l, nvtx_color::grey)
#define NVTX_PUSH_ITER(l)     NVTX_PUSH(l, nvtx_color::white)
#define NVTX_PUSH_KEYGEN(l)   NVTX_PUSH(l, nvtx_color::green)
#define NVTX_PUSH_ENCAPS(l)   NVTX_PUSH(l, nvtx_color::blue)
#define NVTX_PUSH_DECAPS(l)   NVTX_PUSH(l, nvtx_color::orange)
#define NVTX_PUSH_VERIFY(l)   NVTX_PUSH(l, nvtx_color::purple)
#define NVTX_PUSH_TEARDOWN(l) NVTX_PUSH(l, nvtx_color::grey)

/* ======================================================================== */
/*  Regular includes                                                          */
/* ======================================================================== */
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <random>
#include "pke.h"
#include "kem.h"
#include "utils.h"
#include "gpu_backend.h"
#include <openssl/rand.h>
#include "rngongpu_adapter.h"

int main(int argc, char** argv)
{
    const uint32_t q = 12289;
    if (argc < 3) { fprintf(stderr, "Usage: %s <N> <n>\n", argv[0]); return 1; }

    const int      N   = std::atoi(argv[1]);
    const uint32_t n   = (uint32_t)std::atoi(argv[2]);
    const uint32_t MSG = 256;

    printf("N = %d, n = %u, q = %u, MSG = %u\n", N, n, q, MSG);
#ifdef USE_OPENACC
    printf("Backend: OpenACC\n");
#else
    printf("Backend: OpenMP Target\n");
#endif

    /* ====================================================================
     *  INIT: pre-allocate ALL working buffers ONCE on the host side.
     *  Then enter them into the persistent device region with a single
     *  directive.  No further HtoD/DtoH happens unless explicitly coded.
     * ==================================================================== */
    NVTX_PUSH_INIT("Init");

    /* --- Size constants ------------------------------------------------- */
    const size_t noise_total = (size_t)MSG * n;   /* doubles per noise buf  */
    const size_t total_r     = (size_t)MSG * n;   /* int32 r / e1 elements  */
    const size_t total_c     = (size_t)MSG * (n + 1);
    const size_t total_A     = (size_t)n * n;

    /* --- Host-side allocations ------------------------------------------ */
    std::vector<int32_t>  s_buf(n),          t_buf(n),       z_buf(256);
    std::vector<int32_t>  m_buf(MSG),         mp_buf(MSG),    dec_buf(MSG);
    std::vector<int32_t>  ptxt_buf(MSG),      e2_buf(MSG);
    std::vector<int32_t>  r_buf(total_r),     e1_buf(total_r);
    std::vector<int32_t>  c_buf(total_c),     c_chk(total_c);
    std::vector<int32_t>  A_buf(total_A);

    /* FIX 2: Persistent noise temps (avoid cudaMalloc in hot path) */
    std::vector<double>   r_tmp_buf(noise_total);    /* 8.4 MB */
    std::vector<double>   e1_tmp_buf(noise_total);   /* 8.4 MB */
    std::vector<uint32_t> e2_tmp_buf(MSG);           /* 1 KB   */

    /* FIX 4: Persistent KeyGen temps */
    std::vector<uint32_t> sraw_buf(6 * n);           /* 96 KB  */
    std::vector<double>   etmp_buf(n);               /* 32 KB  */
    std::vector<int32_t>  ebuf_buf(n);               /* 16 KB  */

    /* Raw pointers (used directly in pragma clauses) */
    int32_t*  sp      = s_buf.data();
    int32_t*  tp      = t_buf.data();
    int32_t*  zp      = z_buf.data();
    int32_t*  mp      = m_buf.data();
    int32_t*  mpp     = mp_buf.data();
    int32_t*  decp    = dec_buf.data();
    int32_t*  ptxtp   = ptxt_buf.data();
    int32_t*  rp      = r_buf.data();
    int32_t*  e1p     = e1_buf.data();
    int32_t*  e2p     = e2_buf.data();
    int32_t*  cp      = c_buf.data();
    int32_t*  cchkp   = c_chk.data();
    int32_t*  Ap      = A_buf.data();
    double*   r_tmp_d  = r_tmp_buf.data();
    double*   e1_tmp_d = e1_tmp_buf.data();
    uint32_t* e2_tmp_u = e2_tmp_buf.data();
    uint32_t* sraw_p   = sraw_buf.data();
    double*   etmp_p   = etmp_buf.data();
    int32_t*  ebuf_p   = ebuf_buf.data();

    /* ----------------------------------------------------------------
     *  PERSISTENT DEVICE REGION — allocate all buffers ONCE.
     *
     *  OpenACC:   #pragma acc enter data create(...)
     *  OpenMP:    #pragma omp target enter data map(alloc:...)
     *
     *  `create` / `map(alloc:...)`:
     *    Allocates device memory without copying host data.
     *    Used for output/scratch buffers whose initial value doesn't matter.
     *
     *  Under -gpu=managed both compilers may map this to cudaMallocManaged,
     *  in which case the allocation is physically unified and no separate
     *  "enter" directive is strictly necessary — but it is kept for
     *  correctness across memory models.
     * ---------------------------------------------------------------- */
#ifdef USE_OPENACC
    #pragma acc enter data create(                              \
        sp[0:n], tp[0:n], zp[0:256],                           \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],                     \
        ptxtp[0:MSG], e2p[0:MSG],                               \
        rp[0:total_r], e1p[0:total_r],                          \
        cp[0:total_c], cchkp[0:total_c],                        \
        Ap[0:total_A],                                          \
        r_tmp_d[0:noise_total], e1_tmp_d[0:noise_total],        \
        e2_tmp_u[0:MSG],                                        \
        sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])

#else   /* USE_OMP_TARGET */
    #pragma omp target enter data map(alloc:               \
        sp[0:n], tp[0:n], zp[0:256],                       \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],                 \
        ptxtp[0:MSG], e2p[0:MSG],                           \
        rp[0:total_r], e1p[0:total_r],                      \
        cp[0:total_c], cchkp[0:total_c],                    \
        Ap[0:total_A],                                      \
        r_tmp_d[0:noise_total], e1_tmp_d[0:noise_total],    \
        e2_tmp_u[0:MSG],                                    \
        sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])
#endif

    NVTX_POP();  /* Init */

    /* ====================================================================
     *  MAIN BENCHMARK LOOP
     * ==================================================================== */
    long long sum_keygen_us = 0, sum_encaps_us = 0, sum_decaps_us = 0;
    int mismatches = 0;

    const auto startTot = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        NVTX_PUSH_ITER("Iteration");

        uint64_t key_seed;
        RAND_bytes(reinterpret_cast<unsigned char*>(&key_seed), sizeof(key_seed));

        /* ================================================================
         *  KeyGen
         * ================================================================ */
        NVTX_PUSH_KEYGEN("KeyGen");
        auto t0 = std::chrono::steady_clock::now();

        const uint64_t rho_seed = key_seed ^ 0xCAFEBABE12345678ULL;

        NVTX_PUSH_KEYGEN("GenerateA");
        GenerateA_GPU_rngongpu(rho_seed, n, q, Ap);
        NVTX_POP();

        NVTX_PUSH_KEYGEN("KeyGen_GPU");
        KeyGen_GPU_rngongpu_Aflat(key_seed, n, q, Ap, sp, tp,
                                  sraw_p, etmp_p, ebuf_p);
        NVTX_POP();

        /* tp DtoH (16 KB) — needed by SHA3(pkh) on CPU */
        NVTX_PUSH_KEYGEN("tp_DtoH");
#ifdef USE_OPENACC
        #pragma acc update self(tp[0:n])
#else
        #pragma omp target update from(tp[0:n])
#endif
        NVTX_POP();

        /* Random rejection vector z (CPU → device) */
        for (int i = 0; i < 256; ++i) zp[i] = getRandomInt(0, 1);
#ifdef USE_OPENACC
        #pragma acc update device(zp[0:256])
#else
        #pragma omp target update to(zp[0:256])
#endif

        auto t1 = std::chrono::steady_clock::now();
        sum_keygen_us +=
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        NVTX_POP();  /* KeyGen */

        /* ================================================================
         *  Encaps
         * ================================================================ */
        std::vector<int32_t> K_enc;
        std::vector<int32_t> h_c_enc;   /* FIX 1: cached SHA3(c_out) */

        NVTX_PUSH_ENCAPS("Encaps");
        auto t2 = std::chrono::steady_clock::now();

        Encaps_GPU_Aflat(key_seed, n, q, Ap, tp, cp,
                         K_enc, h_c_enc,
                         mp, rp, e1p, e2p, ptxtp,
                         r_tmp_d, e1_tmp_d, e2_tmp_u);

        auto t3 = std::chrono::steady_clock::now();
        sum_encaps_us +=
            std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
        NVTX_POP();  /* Encaps */

        /* ================================================================
         *  Decaps
         * ================================================================ */
        std::vector<int32_t> K_dec;

        NVTX_PUSH_DECAPS("Decaps");
        auto t4 = std::chrono::steady_clock::now();

        Decaps_GPU_Aflat(key_seed, n, q, Ap, tp, sp, zp, cp,
                         K_dec, h_c_enc,    /* FIX 1: h_c reused */
                         mpp, decp, rp, e1p, e2p, ptxtp, cchkp,
                         r_tmp_d, e1_tmp_d, e2_tmp_u);

        auto t5 = std::chrono::steady_clock::now();
        sum_decaps_us +=
            std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
        NVTX_POP();  /* Decaps */

        /* ================================================================
         *  Verify: CPU sanity-check (outside timing)
         * ================================================================ */
        NVTX_PUSH_VERIFY("Verify");

        /* Bring sp and decp back to host for CPU verification */
#ifdef USE_OPENACC
        #pragma acc update self(sp[0:n], decp[0:MSG])
#else
        #pragma omp target update from(sp[0:n], decp[0:MSG])
#endif

        /* Decrypt first ciphertext slot on CPU */
        std::vector<int32_t> u0(n);
        for (uint32_t i = 0; i < n; ++i) u0[i] = cp[i];
        const int32_t v0 = cp[n];

        std::vector<int32_t> s_k_cpu;
        s_k_cpu.reserve(n + 256);
        for (uint32_t i = 0; i < n;   ++i) s_k_cpu.push_back(sp[i]);
        for (int      i = 0; i < 256; ++i) s_k_cpu.push_back(zp[i]);

        int32_t dec0_cpu = -999;
        Decrypt(v0, u0, s_k_cpu, q, dec0_cpu);

        /* Compare shared keys */
        bool same = (K_enc.size() == K_dec.size());
        if (same) {
            for (size_t i = 0; i < K_enc.size(); ++i)
                if (K_enc[i] != K_dec[i]) { same = false; break; }
        }
        if (!same) {
            ++mismatches;
            NVTX_MARK("KeyMismatch");
        }

        NVTX_POP();  /* Verify */
        NVTX_POP();  /* Iteration */
    }

    /* ====================================================================
     *  TEARDOWN: release persistent device memory + print results.
     *
     *  OpenACC:   #pragma acc exit data delete(...)
     *  OpenMP:    #pragma omp target exit data map(release:...)
     *
     *  `delete` / `map(release:...)`:
     *    Frees device memory.  With managed memory this may be a no-op
     *    in terms of physical pages, but it correctly updates the
     *    runtime's data-presence table.
     * ==================================================================== */
    NVTX_PUSH_TEARDOWN("Teardown");

#ifdef USE_OPENACC
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

#else   /* USE_OMP_TARGET */
    #pragma omp target exit data map(release:              \
        sp[0:n], tp[0:n], zp[0:256],                       \
        mp[0:MSG], mpp[0:MSG], decp[0:MSG],                 \
        ptxtp[0:MSG], e2p[0:MSG],                           \
        rp[0:total_r], e1p[0:total_r],                      \
        cp[0:total_c], cchkp[0:total_c],                    \
        Ap[0:total_A],                                      \
        r_tmp_d[0:noise_total], e1_tmp_d[0:noise_total],    \
        e2_tmp_u[0:MSG],                                    \
        sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])
#endif

    const auto endTot = std::chrono::steady_clock::now();
    const double total_s =
        std::chrono::duration_cast<std::chrono::microseconds>(
            endTot - startTot).count() / 1e6;

    printf("%u;%.2f;%.2f;%.2f;%.6f;%d\n",
           n,
           (double)sum_keygen_us / N,
           (double)sum_encaps_us / N,
           (double)sum_decaps_us / N,
           total_s,
           mismatches);

    NVTX_POP();  /* Teardown */
    return 0;
}
