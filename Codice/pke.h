#pragma once
/*  pke.h  –  Public API for LWE PKE GPU kernels  (v4)
 *
 *  Changes vs v3:
 *    OPT-5  uint16_t A  — elements in [0,3328] fit in uint16_t; A_flat,
 *           A_next, DoubleBufA buffers all changed → halves HBM bandwidth.
 *    OPT-6  A^T buffer  — DevicePool::AT holds the row-major transpose of A;
 *           BatchEncrypt re-structured to gang over output element i (4096
 *           CTAs) rather than message j (256 CTAs), reading AT rows once and
 *           sharing them across all 256 j outputs → 256× less A DRAM traffic.
 *    OPT-7  Barrett reduction — barrett_reduce() replaces all % q divisions;
 *           hardware integer division (20–40 cycles) → multiply+shift (~4 cy).
 *    OPT-8  Transposed r layout — r_buf_T[l*MSG+j] / e1_buf_T[l*MSG+j]
 *           instead of [j*n+l]; enables coalesced r reads in restructured
 *           BatchEncrypt (thread j at fixed l accesses consecutive addresses).
 *    OPT-9  int32 accumulation with periodic reduction — inner GEMV loop
 *           uses int32 for 16-iteration tiles, reducing after each tile;
 *           max value per tile = 16 × 3328² = 177M < 2^31.  Halves register
 *           pressure and avoids slow int64 SIMD paths on A100.
 */
#include <cstdint>
#include <vector>

/* ══════════════════════════════════════════════════════════════════════════
 *  Async queue IDs — inline constexpr shared across all TUs via this header.
 * ══════════════════════════════════════════════════════════════════════════ */
inline constexpr int STREAM_A        = 1;
inline constexpr int STREAM_S        = 2;
inline constexpr int STREAM_E        = 3;
inline constexpr int STREAM_PREFETCH = 4;

/* ══════════════════════════════════════════════════════════════════════════
 *  OPT-7: Barrett reduction for the fixed modulus q = 3329
 *
 *  For 64-bit accumulator a in [0, n*(q-1)^2] ≈ [0, 4.54e10] (n=4096):
 *    k=40, m=floor(2^40/q)=330192
 *    r(a) = a − q * floor(a*m / 2^40)
 *    At most one correction needed because floor underestimates by ≤ 2.
 *
 *  For the signed decryption path (a can be slightly negative after subtraction),
 *  the function first normalises a into [0, ∞) before applying Barrett.
 *
 *  Device-callable: #pragma acc routine seq.  The pragma is a no-op on host.
 * ══════════════════════════════════════════════════════════════════════════ */
static constexpr uint64_t BARRETT_M  = 330192ULL;  /* floor(2^40 / 3329) */
static constexpr uint64_t BARRETT_SH = 40ULL;
static constexpr int32_t  Q_VAL      = 3329;

#pragma acc routine seq
inline int32_t barrett_reduce(long long a)
{
    /* Normalise negatives (only occurs in BatchDecrypt subtraction) */
    if (__builtin_expect(a < 0, 0)) {
        long long k = (-a + Q_VAL - 1) / Q_VAL;
        a += k * (long long)Q_VAL;
    }
    uint64_t ua = (uint64_t)a;
    uint64_t q_ = (ua * BARRETT_M) >> BARRETT_SH;
    int32_t  r  = (int32_t)(ua - q_ * (uint64_t)Q_VAL);
    if (r >= Q_VAL) r -= Q_VAL;
    return r;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  DevicePool — persistent scratch-buffer pool (methods in pke_optimized.cpp)
 *
 *  OPT-5: A_next is now uint16_t.
 *  OPT-6: AT (A-transpose, uint16_t, n×n) added.
 *  OPT-8: r_buf_T / e1_buf_T store rounded noise in l-major (transposed) order.
 * ══════════════════════════════════════════════════════════════════════════ */
struct DevicePool {
    /* KeyGen scratch */
    std::vector<uint32_t>  s_raw;       /* 6·n elements  – CBD source words     */
    std::vector<int32_t>   e_buf;       /* n   elements  – rounded Gaussian e   */
    std::vector<double>    e_tmp;       /* n   elements  – raw Gaussian doubles  */

    /* Noise scratch (Encaps + Decaps) */
    std::vector<double>    r_tmp;       /* MSG·n doubles – raw r, j-major        */
    std::vector<double>    e1_tmp;      /* MSG·n doubles – raw e1, j-major       */
    std::vector<int32_t>   r_buf_T;     /* n·MSG int32   – rounded r, l-major   */
    std::vector<int32_t>   e1_buf_T;    /* n·MSG int32   – rounded e1, l-major  */
    std::vector<uint32_t>  e2_tmp;      /* MSG uint32    – raw e2 uniform words  */

    /* OPT-6: A-transpose for coalesced BatchEncrypt */
    std::vector<uint16_t>  AT;          /* n·n uint16_t  – row-major A^T         */

    /* OPT-4 + OPT-5: second A buffer for double-buffered prefetch */
    std::vector<uint16_t>  A_next;      /* n·n uint16_t                          */

    uint32_t n_alloc   = 0;
    uint32_t msg_alloc = 0;

    void init(uint32_t n, uint32_t MSG);
    void release();
    ~DevicePool();
};

/* ── A-matrix generation ─────────────────────────────────────────────────── */
/* OPT-5: A_flat is now uint16_t*                                             */
void GenerateA_GPU_async   (uint64_t rho_seed, uint32_t n, uint32_t q,
                            uint16_t* A_flat, int stream);
void GenerateA_GPU_rngongpu(uint64_t rho_seed, uint32_t n, uint32_t q,
                            uint16_t* A_flat);

/* ── Noise generation ────────────────────────────────────────────────────── */
/* OPT-8: r_buf / e1_buf parameters now receive the l-major (transposed) arrays */
void GPU_GenerateNoise_rngongpu(uint64_t noise_seed, uint32_t n,
                                uint32_t msg_bits,
                                int32_t* r_buf_T, int32_t* e1_buf_T,
                                int32_t* e2_buf,
                                DevicePool& pool);

/* ── KeyGen ──────────────────────────────────────────────────────────────── */
void KeyGen_GPU_rngongpu_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                               uint16_t* A_flat,
                               int32_t* s_out, int32_t* t_out,
                               DevicePool& pool);

/* ── Encrypt ─────────────────────────────────────────────────────────────── */
/* Pre-rounded path — called by Decaps with r_buf_T / e1_buf_T (l-major).
   OPT-5: A_flat is uint16_t. OPT-6: uses pool.AT internally.               */
void BatchEncrypt_GPU_Aflat(uint32_t n, uint32_t q,
                            uint16_t* A_flat, int32_t* t_ptr,
                            int32_t* r_buf_T, int32_t* e1_buf_T,
                            int32_t* e2_buf,
                            int32_t* ptxt_ptr,
                            int32_t* c_out, uint32_t msg_bits,
                            DevicePool& pool);

/* ── Decrypt ─────────────────────────────────────────────────────────────── */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits);
