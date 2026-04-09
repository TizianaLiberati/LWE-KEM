/*  Codice_modificato_v4.cpp  –  LWE-KEM GPU benchmark driver
 *
 *  v4 changes vs v3:
 *    OPT-5: DoubleBufA ping_buf/pong_buf are now std::vector<uint16_t>.
 *    OPT-5: rp/e1p/e2p in main use the transposed pool arrays (r_buf_T /
 *           e1_buf_T) via pool; the local r_buf/e1_buf vectors are removed
 *           because BatchEncrypt_GPU_Aflat now reads from pool.AT and
 *           pool.r_buf_T directly.
 *    All STREAM_* constants come from pke.h — no local redeclaration.
 *
 *  Compile:
 *    nvc++ -std=c++20 -O3 -acc -gpu=cc80,managed,lineinfo -Minfo=accel \
 *          Codice_modificato_v4.cpp pke_optimized.cpp kem.cpp utils.cpp \
 *          -lcurand -lssl -lcrypto -lpthread -o lwe_kem_v4
 *
 *  Run:
 *    ./lwe_kem_v4 <N> <n>
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <random>
#include "pke.h"    /* DevicePool, STREAM_A/S/E/PREFETCH, uint16_t A types */
#include "kem.h"
#include "utils.h"
#include <openssl/rand.h>

/* NOTE: All STREAM_* constants are inline constexpr in pke.h. Do NOT
 *       redeclare them here.                                                */

/* ════════════════════════════════════════════════════════════════════════════
 *  OPT-4 + OPT-5: Double-buffered A manager
 *  Buffer type changed to uint16_t (was int32_t).
 * ════════════════════════════════════════════════════════════════════════════ */
struct DoubleBufA {
    std::vector<uint16_t> ping_buf;
    std::vector<uint16_t> pong_buf;

    uint16_t* ping = nullptr;
    uint16_t* pong = nullptr;

    uint64_t next_rho       = 0;
    uint32_t n_             = 0;
    uint32_t q_             = 0;
    bool     prefetch_valid = false;

    void init(uint32_t n, uint32_t q, uint64_t first_rho_seed) {
        n_ = n;  q_ = q;
        size_t An = (size_t)n * n;
        ping_buf.assign(An, 0);
        pong_buf.assign(An, 0);
        ping = ping_buf.data();
        pong = pong_buf.data();

        #pragma acc enter data create(ping[0:An], pong[0:An])

        GenerateA_GPU_async(first_rho_seed, n, q, ping, STREAM_A);
        #pragma acc wait(STREAM_A)

        next_rho       = first_rho_seed;
        prefetch_valid = false;
    }

    void prefetch_next(uint64_t rho_next) {
        next_rho = rho_next;
        GenerateA_GPU_async(rho_next, n_, q_, pong, STREAM_PREFETCH);
        prefetch_valid = true;
    }

    uint16_t* acquire_next() {
        if (prefetch_valid) {
            #pragma acc wait(STREAM_PREFETCH)
            std::swap(ping, pong);
            prefetch_valid = false;
        } else {
            GenerateA_GPU_async(next_rho, n_, q_, ping, STREAM_A);
            #pragma acc wait(STREAM_A)
        }
        return ping;
    }

    void release() {
        if (!n_) return;
        size_t An = (size_t)n_ * n_;
        #pragma acc wait(STREAM_PREFETCH)
        #pragma acc exit data delete(ping[0:An], pong[0:An])
        n_ = 0;
    }

    ~DoubleBufA() { release(); }
};

/* ════════════════════════════════════════════════════════════════════════════
 *  main
 * ════════════════════════════════════════════════════════════════════════════ */
int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <N> <n>\n";
        return 1;
    }

    const int      N   = std::atoi(argv[1]);
    const uint32_t n   = static_cast<uint32_t>(std::atoi(argv[2]));
    const uint32_t q   = 3329;
    const uint32_t MSG = 256;

    std::cout << "N = " << N << ", n = " << n
              << ", q = " << q << ", MSG = " << MSG << "\n";

    /* Pre-allocate all working buffers ONCE.
       r_buf / e1_buf (n*MSG) removed — these now live in pool.r_buf_T /
       pool.e1_buf_T, which BatchEncrypt_GPU_Aflat reads directly.           */
    const size_t total_c = static_cast<size_t>(MSG) * (n + 1);

    std::vector<int32_t> s_buf   (n,       0);
    std::vector<int32_t> t_buf   (n,       0);
    std::vector<int32_t> z_buf   (256,     0);
    std::vector<int32_t> m_buf   (MSG,     0);
    std::vector<int32_t> mp_buf  (MSG,     0);
    std::vector<int32_t> dec_buf (MSG,     0);
    std::vector<int32_t> ptxt_buf(MSG,     0);
    std::vector<int32_t> e2_buf  (MSG,     0);
    std::vector<int32_t> c_buf   (total_c, 0);
    std::vector<int32_t> c_chk   (total_c, 0);

    int32_t* sp    = s_buf   .data();
    int32_t* tp    = t_buf   .data();
    int32_t* zp    = z_buf   .data();
    int32_t* mp    = m_buf   .data();
    int32_t* mpp   = mp_buf  .data();
    int32_t* decp  = dec_buf .data();
    int32_t* ptxtp = ptxt_buf.data();
    int32_t* e2p   = e2_buf  .data();
    int32_t* cp    = c_buf   .data();
    int32_t* cchkp = c_chk   .data();

    /* Enter working buffers into device memory ONCE */
    #pragma acc enter data create(   \
        sp   [0:n],                  \
        tp   [0:n],                  \
        zp   [0:256],                \
        mp   [0:MSG],                \
        mpp  [0:MSG],                \
        decp [0:MSG],                \
        ptxtp[0:MSG],                \
        e2p  [0:MSG],                \
        cp   [0:total_c],            \
        cchkp[0:total_c])

    /* DevicePool owns all scratch space (r_buf_T, e1_buf_T, AT, etc.) */
    DevicePool pool;
    pool.init(n, MSG);

    /* Convenience aliases to pool's transposed noise buffers */
    int32_t* rp_T  = pool.r_buf_T .data();
    int32_t* e1p_T = pool.e1_buf_T.data();

    /* OPT-4: Fill first A matrix synchronously */
    uint64_t key_seed = 0;
    RAND_bytes(reinterpret_cast<unsigned char*>(&key_seed), sizeof(key_seed));
    uint64_t rho_seed = key_seed ^ 0xCAFEBABE12345678ULL;

    DoubleBufA dba;
    dba.init(n, q, rho_seed);

    long long sum_keygen_us = 0;
    long long sum_encaps_us = 0;
    long long sum_decaps_us = 0;
    int       mismatches    = 0;

    auto startTot = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        RAND_bytes(reinterpret_cast<unsigned char*>(&key_seed), sizeof(key_seed));
        rho_seed = key_seed ^ 0xCAFEBABE12345678ULL;

        /* OPT-4: swap in prefetched A */
        uint16_t* Ap = dba.acquire_next();

        /* ── KeyGen ─────────────────────────────────────────────────────── */
        auto t0 = std::chrono::steady_clock::now();

        /* OPT-6: also builds pool.AT inside KeyGen */
        KeyGen_GPU_rngongpu_Aflat(key_seed, n, q, Ap, sp, tp, pool);

        for (int i = 0; i < 256; ++i) zp[i] = getRandomInt(0, 1);
        #pragma acc update device(zp[0:256])

        auto t1 = std::chrono::steady_clock::now();
        sum_keygen_us +=
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

        /* OPT-4: kick off next-iter A prefetch — overlaps Encaps + Decaps */
        if (k + 1 < N) {
            uint64_t rho_next = 0;
            RAND_bytes(reinterpret_cast<unsigned char*>(&rho_next), sizeof(rho_next));
            dba.prefetch_next(rho_next);
        }

        /* ── Encaps ──────────────────────────────────────────────────────── */
        std::vector<int32_t> K_enc;
        auto t2 = std::chrono::steady_clock::now();

        Encaps_GPU_Aflat(key_seed, n, q, Ap, tp, cp, K_enc,
                         mp, rp_T, e1p_T, e2p, ptxtp,
                         pool);

        auto t3 = std::chrono::steady_clock::now();
        sum_encaps_us +=
            std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();

        /* ── Decaps ──────────────────────────────────────────────────────── */
        std::vector<int32_t> K_dec;
        auto t4 = std::chrono::steady_clock::now();

        Decaps_GPU_Aflat(key_seed, n, q, Ap, tp, sp, zp, cp, K_dec,
                         mpp, decp, rp_T, e1p_T, e2p, ptxtp, cchkp,
                         pool);

        auto t5 = std::chrono::steady_clock::now();
        sum_decaps_us +=
            std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();

        /* ── Correctness check ───────────────────────────────────────────── */
        bool same = (K_enc.size() == K_dec.size());
        if (same)
            for (size_t i = 0; i < K_enc.size(); ++i)
                if (K_enc[i] != K_dec[i]) { same = false; break; }
        if (!same) ++mismatches;
    }

    /* Release in reverse-allocation order */
    pool.release();
    dba.release();

    #pragma acc exit data delete(    \
        sp   [0:n],                  \
        tp   [0:n],                  \
        zp   [0:256],                \
        mp   [0:MSG],                \
        mpp  [0:MSG],                \
        decp [0:MSG],                \
        ptxtp[0:MSG],                \
        e2p  [0:MSG],                \
        cp   [0:total_c],            \
        cchkp[0:total_c])

    auto endTot = std::chrono::steady_clock::now();
    double total_s =
        std::chrono::duration_cast<std::chrono::microseconds>(
            endTot - startTot).count() / 1e6;

    std::cout << "mismatches = " << mismatches << "\n";
    std::cout << n                                      << ";"
              << static_cast<double>(sum_keygen_us) / N << ";"
              << static_cast<double>(sum_encaps_us) / N << ";"
              << static_cast<double>(sum_decaps_us) / N << ";"
              << total_s                                << "\n";

    return 0;
}
