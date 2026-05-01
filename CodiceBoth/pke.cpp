/*  pke.cpp  –  LWE public-key encryption (GPU-accelerated)
 *
 *  Compile-time backend selection:
 *    nvc++ -DUSE_OPENACC    -acc -gpu=cc90,managed ...
 *    nvc++ -DUSE_OMP_TARGET -mp=gpu -gpu=cc90,managed ...
 *
 *  This version avoids NVHPC ICEs like:
 *    "BAD sptr in var_refsym" by not putting array-length variables
 *    in firstprivate clauses on offloaded loops.[web:1][web:3]
 */

#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>
#include "pke.h"
#include "utils.h"
#include "gpu_backend.h"   /* compile-time ACC / OMP-target abstraction */

extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma,
                                     uint64_t seed, uint64_t stream_id);
extern "C" void rngongpu_fill_uniform_u32(uint32_t* d_out, int n,
                                          uint64_t seed, uint64_t stream_id);

/* ===========================================================================
 *  GPU_CiphertextEqual
 *  Returns true when c1[0..len-1] == c2[0..len-1] (device-resident).
 *  A single integer `mismatch` is accumulated with a parallel reduction.
 * =========================================================================== */
bool GPU_CiphertextEqual(const int32_t* c1, const int32_t* c2, size_t len)
{
    int mismatch = 0;

#ifdef USE_OPENACC
    #pragma acc parallel loop gang vector reduction(+:mismatch) \
        present(c1[0:len], c2[0:len])
#else   /* USE_OMP_TARGET */
    #pragma omp target teams distribute parallel for \
        reduction(+:mismatch) \
        is_device_ptr(c1, c2)
#endif
    for (size_t k = 0; k < len; ++k)
        if (c1[k] != c2[k]) mismatch += 1;

    return (mismatch == 0);
}

/* ===========================================================================
 *  GPU_GenerateNoise_rngongpu
 *
 *  Fast path  (r_tmp_d != nullptr):
 *    All three intermediate buffers are pre-allocated and device-resident.
 *    Avoids cudaMalloc/cudaFree inside the critical path.
 *    Two async queues (noise ↔ error) are used when possible.
 *
 *  Slow path  (r_tmp_d == nullptr):
 *    Allocates local std::vectors and creates device copies transiently.
 *    Used by callers that do not maintain a persistent device region.
 * =========================================================================== */
void GPU_GenerateNoise_rngongpu(uint64_t noise_seed, uint32_t n,
                                uint32_t msg_bits,
                                int32_t*  r_buf,    int32_t*  e1_buf,
                                int32_t*  e2_buf,
                                double*   r_tmp_d,  /* nullable */
                                double*   e1_tmp_d, /* nullable */
                                uint32_t* e2_tmp_u) /* nullable */
{
    const size_t total = static_cast<size_t>(msg_bits) * n;
    const size_t mb    = static_cast<size_t>(msg_bits);

    if (r_tmp_d != nullptr)
    {
        /* ----------------------------------------------------------------
         *  FAST PATH – caller owns the persistent device buffers
         * ---------------------------------------------------------------- */

#ifdef USE_OPENACC
        /* Assert residency; expose raw device pointers for the CUDA kernels */
        #pragma acc data present(r_tmp_d[0:total], e1_tmp_d[0:total], \
                                 e2_tmp_u[0:mb],                       \
                                 r_buf[0:total], e1_buf[0:total],      \
                                 e2_buf[0:mb])
        {
            /* Obtain native device pointers for the external CUDA fill calls */
            #pragma acc host_data use_device(r_tmp_d, e1_tmp_d, e2_tmp_u)
            {
                rngongpu_fill_normal(r_tmp_d,  static_cast<int>(total), 1.0,
                                     noise_seed, 0xC0C0);
                rngongpu_fill_normal(e1_tmp_d, static_cast<int>(total), 1.0,
                                     noise_seed, 0xD0D0);
                rngongpu_fill_uniform_u32(e2_tmp_u, static_cast<int>(mb),
                                          noise_seed, 0xE0E0);
            }

            /* Queue 1: round doubles → int32 for r and e1 */
            #pragma acc parallel loop async(1) \
                present(r_tmp_d[0:total], e1_tmp_d[0:total], \
                        r_buf[0:total], e1_buf[0:total])
            for (size_t idx = 0; idx < total; ++idx) {
                double vr = r_tmp_d[idx], ve = e1_tmp_d[idx];
                if (vr >  6.0) vr =  6.0;  if (vr < -6.0) vr = -6.0;
                if (ve >  6.0) ve =  6.0;  if (ve < -6.0) ve = -6.0;
                r_buf[idx]  = static_cast<int32_t>(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
                e1_buf[idx] = static_cast<int32_t>(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
            }

            /* Queue 2: ternary noise for e2 */
            #pragma acc parallel loop async(2) \
                present(e2_tmp_u[0:mb], e2_buf[0:mb])
            for (size_t j = 0; j < mb; ++j)
                e2_buf[j] = static_cast<int32_t>(e2_tmp_u[j] % 3) - 1;

            #pragma acc wait(1, 2)
        }

#else   /* USE_OMP_TARGET */
        /* Expose device pointers for external CUDA kernel calls */
        #pragma omp target data use_device_ptr(r_tmp_d, e1_tmp_d, e2_tmp_u)
        {
            rngongpu_fill_normal(r_tmp_d,  static_cast<int>(total), 1.0,
                                 noise_seed, 0xC0C0);
            rngongpu_fill_normal(e1_tmp_d, static_cast<int>(total), 1.0,
                                 noise_seed, 0xD0D0);
            rngongpu_fill_uniform_u32(e2_tmp_u, static_cast<int>(mb),
                                      noise_seed, 0xE0E0);
        }

        /* Round doubles → int32 (r and e1 combined for memory locality) */
        #pragma omp target teams distribute parallel for \
            is_device_ptr(r_tmp_d, e1_tmp_d, r_buf, e1_buf) nowait
        for (size_t idx = 0; idx < total; ++idx) {
            double vr = r_tmp_d[idx], ve = e1_tmp_d[idx];
            if (vr >  6.0) vr =  6.0;  if (vr < -6.0) vr = -6.0;
            if (ve >  6.0) ve =  6.0;  if (ve < -6.0) ve = -6.0;
            r_buf[idx]  = static_cast<int32_t>(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
            e1_buf[idx] = static_cast<int32_t>(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
        }

        /* Ternary noise e2 – independent, can run concurrently */
        #pragma omp target teams distribute parallel for \
            is_device_ptr(e2_tmp_u, e2_buf) nowait
        for (size_t j = 0; j < mb; ++j)
            e2_buf[j] = static_cast<int32_t>(e2_tmp_u[j] % 3) - 1;

        /* Synchronise both async kernels */
        #pragma omp taskwait

#endif  /* backend */
    }
    else
    {
        /* ----------------------------------------------------------------
         *  SLOW PATH – allocate transient device buffers locally
         * ---------------------------------------------------------------- */
        std::vector<double>   r_tmp(total);
        std::vector<double>   e1_tmp(total);
        std::vector<uint32_t> e2_tmp(mb);
        double*   rp  = r_tmp.data();
        double*   e1p = e1_tmp.data();
        uint32_t* e2p = e2_tmp.data();

#ifdef USE_OPENACC
        #pragma acc data create(rp[0:total], e1p[0:total], e2p[0:mb]) \
                        present(r_buf[0:total], e1_buf[0:total], e2_buf[0:mb])
        {
            #pragma acc host_data use_device(rp, e1p, e2p)
            {
                rngongpu_fill_normal(rp,  static_cast<int>(total), 1.0,
                                     noise_seed, 0xC0C0);
                rngongpu_fill_normal(e1p, static_cast<int>(total), 1.0,
                                     noise_seed, 0xD0D0);
                rngongpu_fill_uniform_u32(e2p, static_cast<int>(mb),
                                          noise_seed, 0xE0E0);
            }

            #pragma acc parallel loop async(1) \
                present(rp[0:total], e1p[0:total], r_buf[0:total], e1_buf[0:total])
            for (size_t idx = 0; idx < total; ++idx) {
                double vr = rp[idx], ve = e1p[idx];
                if (vr >  6.0) vr =  6.0;  if (vr < -6.0) vr = -6.0;
                if (ve >  6.0) ve =  6.0;  if (ve < -6.0) ve = -6.0;
                r_buf[idx]  = static_cast<int32_t>(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
                e1_buf[idx] = static_cast<int32_t>(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
            }

            #pragma acc parallel loop async(2) \
                present(e2p[0:mb], e2_buf[0:mb])
            for (size_t j = 0; j < mb; ++j)
                e2_buf[j] = static_cast<int32_t>(e2p[j] % 3) - 1;

            #pragma acc wait(1, 2)
        }

#else   /* USE_OMP_TARGET */
        /* Use standard OpenMP maps, without vendor-specific 'present' attribute */
        #pragma omp target data map(alloc: rp[0:total], e1p[0:total], e2p[0:mb]) \
                                map(tofrom: r_buf[0:total], e1_buf[0:total], e2_buf[0:mb])
        {
            #pragma omp target data use_device_ptr(rp, e1p, e2p)
            {
                rngongpu_fill_normal(rp,  static_cast<int>(total), 1.0,
                                     noise_seed, 0xC0C0);
                rngongpu_fill_normal(e1p, static_cast<int>(total), 1.0,
                                     noise_seed, 0xD0D0);
                rngongpu_fill_uniform_u32(e2p, static_cast<int>(mb),
                                          noise_seed, 0xE0E0);
            }

            #pragma omp target teams distribute parallel for nowait \
                map(alloc: rp[0:total], e1p[0:total]) \
                map(tofrom: r_buf[0:total], e1_buf[0:total])
            for (size_t idx = 0; idx < total; ++idx) {
                double vr = rp[idx], ve = e1p[idx];
                if (vr >  6.0) vr =  6.0;  if (vr < -6.0) vr = -6.0;
                if (ve >  6.0) ve =  6.0;  if (ve < -6.0) ve = -6.0;
                r_buf[idx]  = static_cast<int32_t>(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
                e1_buf[idx] = static_cast<int32_t>(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
            }

            #pragma omp target teams distribute parallel for nowait \
                map(alloc: e2p[0:mb]) \
                map(tofrom: e2_buf[0:mb])
            for (size_t j = 0; j < mb; ++j)
                e2_buf[j] = static_cast<int32_t>(e2p[j] % 3) - 1;

            #pragma omp taskwait
        }
#endif  /* backend */
    }
}

/* ===========================================================================
 *  GenerateA_GPU_rngongpu
 *
 *  Generates the public matrix A (n×n, mod q) on-device.
 *  araw_ptr is a pure device-side intermediate (never needed on host).
 *  A_flat is written on-device and copied out once so the caller can
 *  keep it warm in the persistent region for all subsequent operations.
 * =========================================================================== */
void GenerateA_GPU_rngongpu(uint64_t rho_seed, uint32_t n, uint32_t q,
                            int32_t* A_flat)
{
    const size_t total = static_cast<size_t>(n) * n;

    std::vector<uint32_t> A_raw(total);
    uint32_t* araw_ptr = A_raw.data();

#ifdef USE_OPENACC
    #pragma acc data create(araw_ptr[0:total]) copyout(A_flat[0:total])
    {
        #pragma acc host_data use_device(araw_ptr)
        {
            rngongpu_fill_uniform_u32(araw_ptr, static_cast<int>(total),
                                      rho_seed, 0xAAAA);
        }

        #pragma acc parallel loop present(araw_ptr[0:total], A_flat[0:total])
        for (size_t idx = 0; idx < total; ++idx)
            A_flat[idx] = static_cast<int32_t>(araw_ptr[idx] % q);
    }

#else   /* USE_OMP_TARGET */
    #pragma omp target data map(alloc: araw_ptr[0:total]) \
                            map(from:  A_flat[0:total])
    {
        #pragma omp target data use_device_ptr(araw_ptr)
        {
            rngongpu_fill_uniform_u32(araw_ptr, static_cast<int>(total),
                                      rho_seed, 0xAAAA);
        }

        #pragma omp target teams distribute parallel for \
            map(alloc: araw_ptr[0:total]) \
            map(from:  A_flat[0:total])
        for (size_t idx = 0; idx < total; ++idx)
            A_flat[idx] = static_cast<int32_t>(araw_ptr[idx] % q);
    }
#endif  /* backend */
}

/* ===========================================================================
 *  KeyGen_GPU_rngongpu_Aflat
 *
 *  Generates secret key s and public key t = A·s + e.
 *  All buffers are caller-provided and device-resident.
 * =========================================================================== */
void KeyGen_GPU_rngongpu_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                               int32_t*  A_flat,
                               int32_t*  s_out,  int32_t* t_out,
                               uint32_t* sraw_p, /* persistent: 6n uint32_t */
                               double*   etmp_p, /* persistent: n double    */
                               int32_t*  ebuf_p) /* persistent: n int32_t   */
{
#ifdef USE_OPENACC
    #pragma acc data present(A_flat[0:static_cast<size_t>(n)*n], \
                             s_out[0:n], t_out[0:n],             \
                             sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])
    {
        /* Phase 1: secret s via sum-of-6 uniform bits → ternary */
        #pragma acc host_data use_device(sraw_p)
        {
            rngongpu_fill_uniform_u32(sraw_p, static_cast<int>(6 * n),
                                      key_seed, 0xA0A0);
        }

        #pragma acc parallel loop present(sraw_p[0:6*n], s_out[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t base = 6 * i;
            int32_t a = static_cast<int32_t>(sraw_p[base+0] % 2)
                      + static_cast<int32_t>(sraw_p[base+1] % 2)
                      + static_cast<int32_t>(sraw_p[base+2] % 2);
            int32_t b = static_cast<int32_t>(sraw_p[base+3] % 2)
                      + static_cast<int32_t>(sraw_p[base+4] % 2)
                      + static_cast<int32_t>(sraw_p[base+5] % 2);
            s_out[i] = a - b;
        }

        /* Phase 2: Gaussian noise → rounded int */
        #pragma acc host_data use_device(etmp_p)
        {
            rngongpu_fill_normal(etmp_p, static_cast<int>(n), 1.0,
                                 key_seed, 0xB0B0);
        }

        #pragma acc parallel loop present(etmp_p[0:n], ebuf_p[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            double ve = etmp_p[i];
            if (ve >  6.0) ve =  6.0;
            if (ve < -6.0) ve = -6.0;
            ebuf_p[i] = static_cast<int32_t>(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
        }

        /* Phase 3: t = A·s + e (mod q) */
        #pragma acc parallel loop gang vector \
            present(A_flat[0:static_cast<size_t>(n)*n], \
                    s_out[0:n], ebuf_p[0:n], t_out[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t j = 0; j < n; ++j)
                acc += static_cast<long long>(
                         A_flat[static_cast<size_t>(i) * n + j])
                     * static_cast<long long>(s_out[j]);
            long long v = (acc + static_cast<long long>(ebuf_p[i]))
                          % static_cast<long long>(q);
            if (v < 0) v += q;
            t_out[i] = static_cast<int32_t>(v);
        }
    }

#else   /* USE_OMP_TARGET */
    /* Phase 1: secret s */
    #pragma omp target data use_device_ptr(sraw_p)
    {
        rngongpu_fill_uniform_u32(sraw_p, static_cast<int>(6 * n),
                                  key_seed, 0xA0A0);
    }

    #pragma omp target teams distribute parallel for \
        is_device_ptr(sraw_p, s_out)
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t base = 6 * i;
        int32_t a = static_cast<int32_t>(sraw_p[base+0] % 2)
                  + static_cast<int32_t>(sraw_p[base+1] % 2)
                  + static_cast<int32_t>(sraw_p[base+2] % 2);
        int32_t b = static_cast<int32_t>(sraw_p[base+3] % 2)
                  + static_cast<int32_t>(sraw_p[base+4] % 2)
                  + static_cast<int32_t>(sraw_p[base+5] % 2);
        s_out[i] = a - b;
    }

    /* Phase 2: Gaussian noise */
    #pragma omp target data use_device_ptr(etmp_p)
    {
        rngongpu_fill_normal(etmp_p, static_cast<int>(n), 1.0,
                             key_seed, 0xB0B0);
    }

    #pragma omp target teams distribute parallel for \
        is_device_ptr(etmp_p, ebuf_p)
    for (uint32_t i = 0; i < n; ++i) {
        double ve = etmp_p[i];
        if (ve >  6.0) ve =  6.0;
        if (ve < -6.0) ve = -6.0;
        ebuf_p[i] = static_cast<int32_t>(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
    }

    /* Phase 3: t = A·s + e (mod q) */
    #pragma omp target teams distribute parallel for \
        is_device_ptr(A_flat, s_out, ebuf_p, t_out)
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += static_cast<long long>(
                     A_flat[static_cast<size_t>(i) * n + j])
                 * static_cast<long long>(s_out[j]);
        long long v = (acc + static_cast<long long>(ebuf_p[i]))
                      % static_cast<long long>(q);
        if (v < 0) v += q;
        t_out[i] = static_cast<int32_t>(v);
    }
#endif  /* backend */
}

/* ===========================================================================
 *  BatchEncrypt_GPU_Aflat
 *
 *  Encrypts msg_bits ciphertexts in parallel.
 *  All input/output buffers must be device-resident (present/is_device_ptr).
 * =========================================================================== */
void BatchEncrypt_GPU_Aflat(uint32_t n, uint32_t q,
                            int32_t* A_flat, int32_t* t_ptr,
                            int32_t* r_buf,  int32_t* e1_buf, int32_t* e2_buf,
                            int32_t* ptxt_ptr,
                            int32_t* c_out, uint32_t msg_bits)
{
    const size_t total_r = static_cast<size_t>(msg_bits) * n;
    const size_t total_c = static_cast<size_t>(msg_bits) * (n + 1);
    const size_t total_A = static_cast<size_t>(n) * n;

#ifdef USE_OPENACC
    #pragma acc parallel loop gang \
        present(A_flat[0:total_A], t_ptr[0:n],               \
                r_buf[0:total_r],  e1_buf[0:total_r],        \
                e2_buf[0:msg_bits], ptxt_ptr[0:msg_bits],    \
                c_out[0:total_c])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        const size_t r_off = static_cast<size_t>(j) * n;
        const size_t c_off = static_cast<size_t>(j) * (n + 1);

        /* u_i = sum_l A[l,i] * r[l] + e1[i]  (mod q) */
        #pragma acc loop vector
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t l = 0; l < n; ++l)
                acc += static_cast<long long>(
                         A_flat[static_cast<size_t>(l) * n + i])
                     * static_cast<long long>(r_buf[r_off + l]);
            long long u = (acc + static_cast<long long>(e1_buf[r_off + i]))
                          % static_cast<long long>(q);
            if (u < 0) u += q;
            c_out[c_off + i] = static_cast<int32_t>(u);
        }

        /* v = t·r + e2 + ptxt  (mod q) */
        long long dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += static_cast<long long>(t_ptr[i])
                 * static_cast<long long>(r_buf[r_off + i]);

        long long vv = (dot
                        + static_cast<long long>(e2_buf[j])
                        + static_cast<long long>(ptxt_ptr[j]))
                       % static_cast<long long>(q);
        if (vv < 0) vv += q;
        c_out[c_off + n] = static_cast<int32_t>(vv);
    }

#else   /* USE_OMP_TARGET */
    #pragma omp target teams distribute \
        is_device_ptr(A_flat, t_ptr, r_buf, e1_buf, e2_buf, ptxt_ptr, c_out)
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        const size_t r_off = static_cast<size_t>(j) * n;
        const size_t c_off = static_cast<size_t>(j) * (n + 1);

        /* u_i = A^T·r + e1 — each thread handles one i */
        #pragma omp parallel for simd
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t l = 0; l < n; ++l)
                acc += static_cast<long long>(
                         A_flat[static_cast<size_t>(l) * n + i])
                     * static_cast<long long>(r_buf[r_off + l]);
            long long u = (acc + static_cast<long long>(e1_buf[r_off + i]))
                          % static_cast<long long>(q);
            if (u < 0) u += q;
            c_out[c_off + i] = static_cast<int32_t>(u);
        }

        /* Dot product with reduction */
        long long dot = 0;
        #pragma omp parallel for simd reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += static_cast<long long>(t_ptr[i])
                 * static_cast<long long>(r_buf[r_off + i]);

        long long vv = (dot
                        + static_cast<long long>(e2_buf[j])
                        + static_cast<long long>(ptxt_ptr[j]))
                       % static_cast<long long>(q);
        if (vv < 0) vv += q;
        c_out[c_off + n] = static_cast<int32_t>(vv);
    }
#endif  /* backend */
}

/* ===========================================================================
 *  BatchDecrypt_GPU
 *
 *  Decrypts msg_bits ciphertexts in parallel.
 * =========================================================================== */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits)
{
    const size_t total_c = static_cast<size_t>(msg_bits) * (n + 1);

#ifdef USE_OPENACC
    #pragma acc parallel loop gang \
        copyin(s_ptr[0:n], c_in[0:total_c]) \
        copyout(decrypt_out[0:msg_bits])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        const size_t c_off = static_cast<size_t>(j) * (n + 1);
        long long dot = 0;

        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += static_cast<long long>(s_ptr[i])
                 * static_cast<long long>(c_in[c_off + i]);

        long long mu = (static_cast<long long>(c_in[c_off + n]) - dot)
                       % static_cast<long long>(q);
        if (mu < 0) mu += q;

        const int32_t bound = static_cast<int32_t>(q) / 4;
        decrypt_out[j] =
            (static_cast<int32_t>(mu) <= bound ||
             static_cast<int32_t>(mu) >= static_cast<int32_t>(q) - bound)
            ? 0 : static_cast<int32_t>(q / 2);
    }

#else   /* USE_OMP_TARGET */
    #pragma omp target teams distribute parallel for \
        map(to:   s_ptr[0:n], c_in[0:total_c]) \
        map(from: decrypt_out[0:msg_bits])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        const size_t c_off = static_cast<size_t>(j) * (n + 1);
        long long dot = 0;

        #pragma omp simd reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += static_cast<long long>(s_ptr[i])
                 * static_cast<long long>(c_in[c_off + i]);

        long long mu = (static_cast<long long>(c_in[c_off + n]) - dot)
                       % static_cast<long long>(q);
        if (mu < 0) mu += q;

        const int32_t bound = static_cast<int32_t>(q) / 4;
        decrypt_out[j] =
            (static_cast<int32_t>(mu) <= bound ||
             static_cast<int32_t>(mu) >= static_cast<int32_t>(q) - bound)
            ? 0 : static_cast<int32_t>(q / 2);
    }
#endif  /* backend */
}

/* ===========================================================================
 *  Decrypt  (CPU-facing single-ciphertext helper)
 *
 *  Used by the CPU verification path in the benchmark driver.
 * =========================================================================== */
void Decrypt(int32_t v_i, const std::vector<int32_t>& u,
             const std::vector<int32_t>& s_k, uint32_t q, int32_t& decrypt_i)
{
    const size_t n = u.size();
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n);
    long long dot = 0;
    const int32_t* up = u.data();
    const int32_t* sp = s.data();

#ifdef USE_OPENACC
    #pragma acc parallel loop reduction(+:dot) copyin(up[0:n], sp[0:n])
    for (size_t i = 0; i < n; ++i)
        dot += static_cast<long long>(sp[i]) * static_cast<long long>(up[i]);

#else   /* USE_OMP_TARGET */
    #pragma omp target teams distribute parallel for \
        reduction(+:dot) \
        map(to: up[0:n], sp[0:n])
    for (size_t i = 0; i < n; ++i)
        dot += static_cast<long long>(sp[i]) * static_cast<long long>(up[i]);
#endif

    long long r = (static_cast<long long>(v_i) - dot) % static_cast<long long>(q);
    if (r < 0) r += q;
    const int32_t mu    = static_cast<int32_t>(r);
    const int32_t bound = static_cast<int32_t>(q) / 4;
    decrypt_i = (mu <= bound || mu >= static_cast<int32_t>(q) - bound)
                  ? 0 : static_cast<int32_t>(q / 2);
}
