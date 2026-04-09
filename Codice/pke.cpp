/*  pke_optimized.cpp  –  LWE PKE GPU kernels  (v4 – nsys-guided pass)
 *
 *  NEW IN v4 (on top of v3 OPT-1 .. OPT-4):
 *
 *  OPT-5  uint16_t A storage
 *    A_flat, A_next, DoubleBufA buffers changed from int32_t to uint16_t.
 *    Elements are in [0, q-1] = [0, 3328], which fits in uint16_t (max 65535).
 *    Halves A memory: 64 MB → 32 MB at n=4096.
 *    All GEMV/GEMM kernels cast uint16_t → int32_t on the fly (zero extra cost).
 *
 *  OPT-6  A^T buffer + restructured BatchEncrypt
 *    pool.AT[i*n+l] = A[l*n+i] (row-major transpose, uint16_t).
 *    AT is built once per KeyGen call on STREAM_A alongside the regular A fill.
 *    BatchEncrypt_GPU_Aflat now uses:
 *      • gang over i = 0..n-1  (4096 CTAs; each computes one output row of c)
 *      • vector over j = 0..MSG-1  (256 threads per CTA)
 *      • sequential over l = 0..n-1  (inner dot product)
 *    Thread j reads AT[i*n+l] (broadcast within CTA, same address) and
 *    r_buf_T[l*MSG+j] (stride-1 coalesced for varying j).
 *    A^T rows are loaded once per gang → 4096 × 32 KB = 32 MB total vs 16 GB before.
 *
 *  OPT-7  Barrett reduction for q = 3329
 *    barrett_reduce() (declared inline in pke.h) replaces all % (long long)q.
 *    Hardware integer division: ~20-40 cycles.  Barrett: ~4 cycles (mul+shift).
 *    Applied in KeyGen GEMV, BatchEncrypt, BatchDecrypt inner loops.
 *
 *  OPT-8  Transposed r/e1 layout
 *    GPU_GenerateNoise_rngongpu fills r_buf_T[l*MSG+j] and e1_buf_T[l*MSG+j]
 *    (l-major) instead of [j*n+l] (j-major).
 *    Enables stride-1 access for thread j at fixed l in the restructured encrypt.
 *    No extra transpose kernel needed — the rounding pass writes transposed directly.
 *
 *  OPT-9  int32 accumulation with periodic modular reduction
 *    GEMV inner loop splits into 16-element tiles and reduces mod q after each tile.
 *    Max tile accumulator = 16 × 3328^2 = 177,209,344 < 2^31.
 *    Avoids int64 SIMD paths on A100 (int64 MACs are ~2× slower than int32).
 *
 *  Compile:
 *    nvc++ -std=c++20 -O3 -acc -gpu=cc80,managed,lineinfo -Minfo=accel \
 *          -o lwe_kem pke_optimized.cpp kem.cpp utils.cpp \
 *          -lcurand -lssl -lcrypto
 */

#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <cassert>
#include "pke.h"
#include "utils.h"
#include "xorshift.h"

/* STREAM_* constants come from pke.h — do NOT redeclare here. */

extern "C" void rngongpu_fill_normal      (double*   d_out, int n, double sigma,
                                           uint64_t seed, uint64_t stream_id);
extern "C" void rngongpu_fill_uniform_u32 (uint32_t* d_out, int n,
                                           uint64_t seed, uint64_t stream_id);

/* ══════════════════════════════════════════════════════════════════════════
 *  DevicePool — out-of-line method definitions
 * ══════════════════════════════════════════════════════════════════════════ */
void DevicePool::init(uint32_t n, uint32_t MSG)
{
    n_alloc   = n;
    msg_alloc = MSG;

    /* Round to 256-byte (64 × int32 or 128 × uint16) cache-line boundaries. */
    auto pad16 = [](size_t e){ return (e + 127) & ~size_t(127); }; /* uint16 */
    auto pad32 = [](size_t e){ return (e + 63)  & ~size_t(63);  }; /* int32  */
    auto pad64 = [](size_t e){ return (e + 31)  & ~size_t(31);  }; /* double */

    size_t sn   = pad32(6 * n);
    size_t en   = pad32(n);
    size_t edn  = pad64(n);
    size_t rdn  = pad64((size_t)MSG * n);   /* r_tmp / e1_tmp (j-major double) */
    size_t rtn  = pad32((size_t)n * MSG);   /* r_buf_T / e1_buf_T (l-major)    */
    size_t e2n  = pad32(MSG);
    size_t An   = pad16((size_t)n * n);     /* AT and A_next (uint16_t)         */

    s_raw   .assign(sn,  0);
    e_buf   .assign(en,  0);
    e_tmp   .assign(edn, 0.0);
    r_tmp   .assign(rdn, 0.0);
    e1_tmp  .assign(rdn, 0.0);
    r_buf_T .assign(rtn, 0);
    e1_buf_T.assign(rtn, 0);
    e2_tmp  .assign(e2n, 0);
    AT      .assign(An,  0);
    A_next  .assign(An,  0);

    uint32_t*  sp   = s_raw   .data();
    int32_t*   ep   = e_buf   .data();
    double*    etp  = e_tmp   .data();
    double*    rp   = r_tmp   .data();
    double*    e1p  = e1_tmp  .data();
    int32_t*   rtp  = r_buf_T .data();
    int32_t*   e1tp = e1_buf_T.data();
    uint32_t*  e2p  = e2_tmp  .data();
    uint16_t*  ATp  = AT      .data();
    uint16_t*  ANp  = A_next  .data();

    #pragma acc enter data create(                     \
        sp  [0:sn],   ep  [0:en],   etp [0:edn],      \
        rp  [0:rdn],  e1p [0:rdn],                    \
        rtp [0:rtn],  e1tp[0:rtn],  e2p [0:e2n],      \
        ATp [0:An],   ANp [0:An])
}

void DevicePool::release()
{
    if (!n_alloc) return;

    auto pad16 = [](size_t e){ return (e + 127) & ~size_t(127); };
    auto pad32 = [](size_t e){ return (e + 63)  & ~size_t(63);  };
    auto pad64 = [](size_t e){ return (e + 31)  & ~size_t(31);  };

    size_t sn  = pad32(6 * n_alloc);
    size_t en  = pad32(n_alloc);
    size_t edn = pad64(n_alloc);
    size_t rdn = pad64((size_t)msg_alloc * n_alloc);
    size_t rtn = pad32((size_t)n_alloc * msg_alloc);
    size_t e2n = pad32(msg_alloc);
    size_t An  = pad16((size_t)n_alloc * n_alloc);

    uint32_t*  sp   = s_raw   .data();
    int32_t*   ep   = e_buf   .data();
    double*    etp  = e_tmp   .data();
    double*    rp   = r_tmp   .data();
    double*    e1p  = e1_tmp  .data();
    int32_t*   rtp  = r_buf_T .data();
    int32_t*   e1tp = e1_buf_T.data();
    uint32_t*  e2p  = e2_tmp  .data();
    uint16_t*  ATp  = AT      .data();
    uint16_t*  ANp  = A_next  .data();

    #pragma acc exit data delete(                      \
        sp  [0:sn],   ep  [0:en],   etp [0:edn],      \
        rp  [0:rdn],  e1p [0:rdn],                    \
        rtp [0:rtn],  e1tp[0:rtn],  e2p [0:e2n],      \
        ATp [0:An],   ANp [0:An])

    n_alloc = 0;
}

DevicePool::~DevicePool() { release(); }

/* ──────────────────────────────────────────────────────────────────────────
 *  Device routine: deterministic A element (xorshift, device-callable)
 *  OPT-5: returns uint16_t (was int32_t)
 * ────────────────────────────────────────────────────────────────────────── */
#pragma acc routine seq
static inline uint16_t A_elem(uint64_t key_seed, int i, int j, int n, int q)
{
    uint64_t st = key_seed
                ^ ((uint64_t)i * 0x9E3779B97F4A7C15ULL)
                ^ ((uint64_t)j * 0x517CC1B727220A95ULL);
    return (uint16_t)xorshift_uniform(&st, q - 1);
}

/* ══════════════════════════════════════════════════════════════════════════
 *  GenerateA_GPU_async
 *
 *  OPT-5: A_flat is now uint16_t* (halved memory traffic).
 *  OPT-6: Also fills pool.AT (the transpose) on the same stream.
 *          AT[i*n+j] = A[j*n+i], built with a second collapse(2) loop
 *          that runs immediately after A is written — both loops share the
 *          same async queue so no extra synchronisation is needed.
 * ══════════════════════════════════════════════════════════════════════════ */
void GenerateA_GPU_async(uint64_t rho_seed, uint32_t n, uint32_t q,
                         uint16_t* A_flat, int stream)
{
    size_t total = (size_t)n * n;

    /* Pass A_flat to the transpose kernel later — get AT from caller's pool
       only if a pool ptr is available.  For the simple async wrapper (used by
       DoubleBufA which doesn't carry a pool), we skip the transpose; the main
       KeyGen path calls the pool-aware overload below.                       */
    #pragma acc parallel loop gang vector collapse(2) async(stream) \
        firstprivate(rho_seed, n, q) present(A_flat[0:total])
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            A_flat[(size_t)i * n + j] = A_elem(rho_seed, (int)i, (int)j, (int)n, (int)q);
}

/* Pool-aware overload: fills A_flat AND pool.AT on the same stream. */
static void GenerateA_and_AT_async(uint64_t rho_seed, uint32_t n, uint32_t q,
                                   uint16_t* A_flat, uint16_t* AT_flat, int stream)
{
    size_t total = (size_t)n * n;

    /* Fill A */
    #pragma acc parallel loop gang vector collapse(2) async(stream) \
        firstprivate(rho_seed, n, q) present(A_flat[0:total])
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            A_flat[(size_t)i * n + j] = A_elem(rho_seed, (int)i, (int)j, (int)n, (int)q);

    /* Fill A^T — same queue, runs after A is complete */
    #pragma acc parallel loop gang vector collapse(2) async(stream) \
        firstprivate(n) present(A_flat[0:total], AT_flat[0:total])
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            AT_flat[(size_t)i * n + j] = A_flat[(size_t)j * n + i];
}

/* Convenience synchronous wrapper */
void GenerateA_GPU_rngongpu(uint64_t rho_seed, uint32_t n, uint32_t q,
                            uint16_t* A_flat)
{
    GenerateA_GPU_async(rho_seed, n, q, A_flat, STREAM_A);
    #pragma acc wait(STREAM_A)
}

/* ══════════════════════════════════════════════════════════════════════════
 *  GPU_GenerateNoise_rngongpu
 *
 *  OPT-8: Output is now r_buf_T[l*MSG+j] and e1_buf_T[l*MSG+j] (l-major).
 *  The rounding kernel writes the transposed layout directly — no extra
 *  transpose kernel is needed.
 * ══════════════════════════════════════════════════════════════════════════ */
void GPU_GenerateNoise_rngongpu(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                                int32_t* r_buf_T, int32_t* e1_buf_T,
                                int32_t* e2_buf,
                                DevicePool& pool)
{
    size_t total = (size_t)msg_bits * n;   /* MSG × n elements in j-major double */
    size_t rtn   = (size_t)n * msg_bits;   /* n × MSG elements in l-major int32  */

    auto pad64 = [](size_t e){ return (e + 31) & ~size_t(31); };
    auto pad32 = [](size_t e){ return (e + 63) & ~size_t(63); };
    size_t rdn  = pad64(total);
    size_t rtn2 = pad32(rtn);
    size_t e2n  = ((size_t)msg_bits + 63) & ~size_t(63);

    double*   rp_tmp  = pool.r_tmp .data();
    double*   e1p_tmp = pool.e1_tmp.data();
    uint32_t* e2p_tmp = pool.e2_tmp.data();

    /* Fill raw doubles on device via RNGonGPU */
    #pragma acc host_data use_device(rp_tmp, e1p_tmp, e2p_tmp)
    {
        rngongpu_fill_normal      (rp_tmp,  (int)total,    1.0, noise_seed, 0xC0C0);
        rngongpu_fill_normal      (e1p_tmp, (int)total,    1.0, noise_seed, 0xD0D0);
        rngongpu_fill_uniform_u32 (e2p_tmp, (int)msg_bits, noise_seed,      0xE0E0);
    }

    /* OPT-8: Round double→int32 AND write in transposed l-major layout.
       r_tmp[j*n+l] (j-major) → r_buf_T[l*MSG+j] (l-major).
       Two async queues: one for r+e1, one for e2.                           */
    #pragma acc parallel loop gang vector collapse(2) async(1) \
        firstprivate(n, msg_bits) \
        present(rp_tmp[0:rdn], e1p_tmp[0:rdn], r_buf_T[0:rtn2], e1_buf_T[0:rtn2])
    for (uint32_t l = 0; l < n; ++l) {
        for (uint32_t j = 0; j < msg_bits; ++j) {
            double vr = rp_tmp [(size_t)j * n + l];
            double ve = e1p_tmp[(size_t)j * n + l];
            if (vr >  6.0) vr =  6.0; if (vr < -6.0) vr = -6.0;
            if (ve >  6.0) ve =  6.0; if (ve < -6.0) ve = -6.0;
            /* Write l-major: r_buf_T[l*MSG+j] */
            r_buf_T [(size_t)l * msg_bits + j] = (int32_t)(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
            e1_buf_T[(size_t)l * msg_bits + j] = (int32_t)(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
        }
    }

    #pragma acc parallel loop async(2) \
        present(e2p_tmp[0:e2n], e2_buf[0:msg_bits])
    for (uint32_t j = 0; j < msg_bits; ++j)
        e2_buf[j] = (int32_t)(e2p_tmp[j] % 7) - 3;

    #pragma acc wait(1, 2)
}

/* ══════════════════════════════════════════════════════════════════════════
 *  KeyGen_GPU_rngongpu_Aflat
 *
 *  OPT-5: A_flat is uint16_t.
 *  OPT-6: Builds pool.AT alongside A on STREAM_A.
 *  OPT-7: barrett_reduce() replaces % q in GEMV.
 *  OPT-9: int32 accumulation with 16-element tile reduction.
 * ══════════════════════════════════════════════════════════════════════════ */
void KeyGen_GPU_rngongpu_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                               uint16_t* A_flat,
                               int32_t* s_out, int32_t* t_out,
                               DevicePool& pool)
{
    auto pad16 = [](size_t e){ return (e + 127) & ~size_t(127); };
    auto pad32 = [](size_t e){ return (e + 63)  & ~size_t(63);  };
    auto pad64 = [](size_t e){ return (e + 31)  & ~size_t(31);  };

    size_t sn  = pad32(6 * n);
    size_t en  = pad32(n);
    size_t edn = pad64(n);
    size_t An  = pad16((size_t)n * n);

    uint32_t* sraw_ptr = pool.s_raw.data();
    double*   ep_tmp   = pool.e_tmp.data();
    int32_t*  ebuf_ptr = pool.e_buf.data();
    uint16_t* ATp      = pool.AT.data();

    /* OPT-6: build A^T concurrently with A generation on STREAM_A.
       The transpose is on the same queue so it runs after A is ready,
       with no additional barrier.                                            */
    GenerateA_and_AT_async(key_seed, n, q, A_flat, ATp, STREAM_A);

    /* OPT-3: s_raw and e_tmp fill on independent streams */
    #pragma acc host_data use_device(sraw_ptr, ep_tmp)
    {
        rngongpu_fill_uniform_u32(sraw_ptr, (int)(6 * n), key_seed, 0xA0A0);
        rngongpu_fill_normal     (ep_tmp,   (int)n, 2.3,  key_seed, 0xB0B0);
    }

    /* Compute s_out on STREAM_S */
    #pragma acc parallel loop async(STREAM_S) \
        present(sraw_ptr[0:sn], s_out[0:n])
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t base = 6 * i;
        int32_t a = (int32_t)(sraw_ptr[base+0] % 2)
                  + (int32_t)(sraw_ptr[base+1] % 2)
                  + (int32_t)(sraw_ptr[base+2] % 2);
        int32_t b = (int32_t)(sraw_ptr[base+3] % 2)
                  + (int32_t)(sraw_ptr[base+4] % 2)
                  + (int32_t)(sraw_ptr[base+5] % 2);
        s_out[i] = a - b;
    }

    /* Round e_tmp → ebuf on STREAM_E */
    #pragma acc parallel loop async(STREAM_E) \
        present(ep_tmp[0:edn], ebuf_ptr[0:en])
    for (uint32_t i = 0; i < n; ++i) {
        double ve  = ep_tmp[i];
        double bnd = 6.0 * 2.3;
        if (ve >  bnd) ve =  bnd;
        if (ve < -bnd) ve = -bnd;
        ebuf_ptr[i] = (int32_t)(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
    }

    /* Barrier: wait for A+AT (STREAM_A), s (STREAM_S), e (STREAM_E) */
    #pragma acc wait(STREAM_A, STREAM_S, STREAM_E)

    /* OPT-7 + OPT-9: GEMV with int32 tiles + Barrett reduction */
    #pragma acc parallel loop gang vector firstprivate(n) \
        present(A_flat[0:An], s_out[0:n], ebuf_ptr[0:en], t_out[0:n])
    for (uint32_t i = 0; i < n; ++i) {
        int32_t acc32 = 0;
        /* OPT-9: process 16-element tiles in int32, reduce after each tile */
        for (uint32_t l = 0; l < n; l += 16) {
            uint32_t lend = (l + 16 < n) ? l + 16 : n;
            for (uint32_t j = l; j < lend; ++j)
                acc32 += (int32_t)A_flat[(size_t)i * n + j] * s_out[j];
            /* Keep accumulator in [-q+1, q*(q-1)*16] range */
            acc32 = barrett_reduce((long long)acc32);
        }
        /* OPT-7: final reduction fusing e addition */
        t_out[i] = barrett_reduce((long long)acc32 + (long long)ebuf_ptr[i]);
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 *  BatchEncrypt_GPU_Aflat
 *
 *  OPT-5: A_flat is uint16_t.
 *  OPT-6: Restructured: gang over i (n=4096), vector over j (MSG=256).
 *         Reads pool.AT[i*n+l] (broadcast within CTA from coalesced AT rows)
 *         and r_buf_T[l*MSG+j] (stride-1 coalesced for varying j).
 *         A^T rows are loaded exactly once per gang → 32 MB total vs 16 GB.
 *  OPT-7: barrett_reduce() replaces % q.
 *  OPT-9: int32 tile accumulation.
 * ══════════════════════════════════════════════════════════════════════════ */
void BatchEncrypt_GPU_Aflat(uint32_t n, uint32_t q,
                            uint16_t* A_flat, int32_t* t_ptr,
                            int32_t* r_buf_T, int32_t* e1_buf_T,
                            int32_t* e2_buf,
                            int32_t* ptxt_ptr,
                            int32_t* c_out, uint32_t msg_bits,
                            DevicePool& pool)
{
    auto pad16 = [](size_t e){ return (e + 127) & ~size_t(127); };
    auto pad32 = [](size_t e){ return (e + 63)  & ~size_t(63);  };

    size_t An    = pad16((size_t)n * n);
    size_t rtn   = pad32((size_t)n * msg_bits);
    size_t total_c = (size_t)msg_bits * (n + 1);

    uint16_t* ATp = pool.AT.data();

    /* OPT-6: gang over i (output element), vector over j (message index).
       Each gang i:
         - reads AT[i*n+0..n-1] once (coalesced row read, shared across j)
         - for each j: computes u[j] = AT-row·r_T[:,j] + e1_T[l,j] (mod q)
         - writes c[j*(n+1)+i] for all j
         - computes v[j] = t·r_T[:,j] + e2[j] + ptxt[j] (mod q) → c[j*(n+1)+n]  */
    #pragma acc parallel loop gang firstprivate(n, msg_bits) \
        present(ATp[0:An], t_ptr[0:n],                       \
                r_buf_T[0:rtn], e1_buf_T[0:rtn],             \
                e2_buf[0:msg_bits], ptxt_ptr[0:msg_bits],    \
                c_out[0:total_c])
    for (uint32_t i = 0; i < n; ++i)
    {
        /* u[j] = (AT[i] · r_T[:,j] + e1_T[i,j]) mod q  */
        #pragma acc loop vector
        for (uint32_t j = 0; j < msg_bits; ++j) {
            int32_t acc = 0;
            /* OPT-9: 16-element tiles in int32 */
            for (uint32_t l = 0; l < n; l += 16) {
                uint32_t lend = (l + 16 < n) ? l + 16 : n;
                for (uint32_t ll = l; ll < lend; ++ll)
                    /* AT[i*n+ll]: broadcast (same for all j); r_T[ll*MSG+j]: coalesced */
                    acc += (int32_t)ATp[(size_t)i * n + ll]
                         * r_buf_T[(size_t)ll * msg_bits + j];
                acc = barrett_reduce((long long)acc);
            }
            int32_t u = barrett_reduce((long long)acc + (long long)e1_buf_T[(size_t)i * msg_bits + j]);
            c_out[(size_t)j * (n + 1) + i] = u;
        }
    }

    /* v[j] = (t · r_T[:,j] + e2[j] + ptxt[j]) mod q   — gang over j */
    #pragma acc parallel loop gang vector firstprivate(n, msg_bits) \
        present(t_ptr[0:n], r_buf_T[0:rtn],                        \
                e2_buf[0:msg_bits], ptxt_ptr[0:msg_bits],          \
                c_out[0:total_c])
    for (uint32_t j = 0; j < msg_bits; ++j) {
        int32_t dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t l = 0; l < n; ++l)
            dot += t_ptr[l] * r_buf_T[(size_t)l * msg_bits + j];
        c_out[(size_t)j * (n + 1) + n] =
            barrett_reduce((long long)dot + e2_buf[j] + ptxt_ptr[j]);
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 *  BatchDecrypt_GPU
 *
 *  OPT-7: barrett_reduce() replaces % q.
 * ══════════════════════════════════════════════════════════════════════════ */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits)
{
    size_t total_c = (size_t)msg_bits * (n + 1);

    #pragma acc parallel loop gang firstprivate(n, msg_bits) \
        present(s_ptr[0:n], c_in[0:total_c], decrypt_out[0:msg_bits])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t c_off = (size_t)j * (n + 1);
        int32_t dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += s_ptr[i] * c_in[c_off + i];

        /* OPT-7: Barrett for the subtraction path (can be negative) */
        int32_t mu = barrett_reduce((long long)c_in[c_off + n] - (long long)dot);

        int32_t bound = Q_VAL / 4;
        decrypt_out[j] = (mu <= bound || mu >= Q_VAL - bound) ? 0 : (Q_VAL / 2);
    }
}
