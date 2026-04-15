/*  pke.cpp  –  LWE public-key encryption (GPU-accelerated)
 *
 *  OPTIMIZATIONS in this version (report13 → report14):
 *
 *  ┌──────────────────────────────────────────────────────────────────────┐
 *  │  FIX 2: PERSISTENT TEMP BUFFERS IN GPU_GenerateNoise_rngongpu        │
 *  │                                                                      │
 *  │  Previous version (report13):                                        │
 *  │    Each call to GPU_GenerateNoise_rngongpu allocated on the heap:    │
 *  │      std::vector<double>   r_tmp(total)   = 256*4096*8 = 8.4 MB    │
 *  │      std::vector<double>   e1_tmp(total)  = 8.4 MB                  │
 *  │      std::vector<uint32_t> e2_tmp(msg_bits)= 1 KB                   │
 *  │    These triggered cudaMalloc + cudaFree on EVERY call.              │
 *  │    Called 2× per iter (Encaps+Decaps) = 4 heavy alloc/free cycles.  │
 *  │    The trace shows cudaFree inside the Encaps timing window.         │
 *  │                                                                      │
 *  │  FIX: Accept pre-allocated device pointers as parameters.            │
 *  │    When non-null: use the persistent buffers (zero allocation cost). │
 *  │    When null:     fall back to local allocation (for Encaps_GPU).    │
 *  │                                                                      │
 *  │  The caller (Codice_modificato.cpp) adds these to the persistent     │
 *  │  acc enter data region alongside all other working buffers.          │
 *  └──────────────────────────────────────────────────────────────────────┘
 *
 *  ┌──────────────────────────────────────────────────────────────────────┐
 *  │  FIX 4: PERSISTENT TEMP BUFFERS IN KeyGen_GPU_rngongpu_Aflat         │
 *  │                                                                      │
 *  │  Same issue: s_raw(6n uint32_t), e_tmp(n double), e_buf(n int32_t)  │
 *  │  were allocated per call. Now accepted as parameters from the        │
 *  │  persistent region.                                                  │
 *  └──────────────────────────────────────────────────────────────────────┘
 */
#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>
#include "pke.h"
#include "utils.h"
#include "xorshift.h"

extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma,
                                     uint64_t seed, uint64_t stream_id);
extern "C" void rngongpu_fill_uniform_u32(uint32_t* d_out, int n,
                                          uint64_t seed, uint64_t stream_id);

#pragma acc routine seq
static inline int32_t A_elem(uint64_t key_seed, int i, int j, int n, int q)
{
    uint64_t st = key_seed
                ^ ((uint64_t)i * 0x9E3779B97F4A7C15ULL)
                ^ ((uint64_t)j * 0x517CC1B727220A95ULL);
    return xorshift_uniform(&st, q - 1);
}

/* ================================================================== */
/*  GPU_CiphertextEqual — returns bool (4 bytes DtoH)                  */
/* ================================================================== */
bool GPU_CiphertextEqual(const int32_t* c1, const int32_t* c2, size_t len)
{
    int mismatch = 0;
    #pragma acc parallel loop gang vector reduction(+:mismatch) \
        present(c1[0:len], c2[0:len]) firstprivate(len)
    for (size_t k = 0; k < len; ++k)
        if (c1[k] != c2[k]) mismatch += 1;
    return (mismatch == 0);
}

/* ================================================================== */
/*  GPU_GenerateNoise_rngongpu                                          */
/*                                                                     */
/*  FIX 2: accepts pre-allocated device pointers for the double/u32   */
/*  intermediates. When non-null these are present on device; the      */
/*  local std::vector fallback is only used for the Encaps_GPU (on-   */
/*  the-fly A) variant that does not participate in the persistent     */
/*  device region.                                                     */
/* ================================================================== */
void GPU_GenerateNoise_rngongpu(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                                int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                                double*   r_tmp_d,   /* nullable — persistent device ptr */
                                double*   e1_tmp_d,  /* nullable */
                                uint32_t* e2_tmp_u)  /* nullable */
{
    size_t total = (size_t)msg_bits * n;

    if (r_tmp_d != nullptr) {
        /* ---- FAST PATH: use caller-provided persistent device buffers ---- */
        #pragma acc data present(r_tmp_d[0:total], e1_tmp_d[0:total], \
                                 e2_tmp_u[0:msg_bits],                 \
                                 r_buf[0:total], e1_buf[0:total],       \
                                 e2_buf[0:msg_bits])
        {
            #pragma acc host_data use_device(r_tmp_d, e1_tmp_d, e2_tmp_u)
            {
                rngongpu_fill_normal(r_tmp_d,  (int)total,    1.0, noise_seed, 0xC0C0);
                rngongpu_fill_normal(e1_tmp_d, (int)total,    1.0, noise_seed, 0xD0D0);
                rngongpu_fill_uniform_u32(e2_tmp_u, (int)msg_bits, noise_seed, 0xE0E0);
            }

            #pragma acc parallel loop async(1) \
                present(r_tmp_d[0:total], e1_tmp_d[0:total], \
                        r_buf[0:total], e1_buf[0:total])
            for (size_t idx = 0; idx < total; ++idx) {
                double vr = r_tmp_d[idx],  ve = e1_tmp_d[idx];
                if (vr >  6.0) vr =  6.0; if (vr < -6.0) vr = -6.0;
                if (ve >  6.0) ve =  6.0; if (ve < -6.0) ve = -6.0;
                r_buf[idx]  = (int32_t)(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
                e1_buf[idx] = (int32_t)(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
            }

            #pragma acc parallel loop async(2) \
                present(e2_tmp_u[0:msg_bits], e2_buf[0:msg_bits])
            for (uint32_t j = 0; j < msg_bits; ++j)
                e2_buf[j] = (int32_t)(e2_tmp_u[j] % 3) - 1;

            #pragma acc wait(1,2)
        }
    } else {
        /* ---- SLOW PATH: allocate locally (fallback for non-persistent callers) ---- */
        std::vector<double>   r_tmp(total);
        std::vector<double>   e1_tmp(total);
        std::vector<uint32_t> e2_tmp(msg_bits);
        double*   rp  = r_tmp.data();
        double*   e1p = e1_tmp.data();
        uint32_t* e2p = e2_tmp.data();

        #pragma acc data create(rp[0:total], e1p[0:total], e2p[0:msg_bits]) \
                        present(r_buf[0:total], e1_buf[0:total], e2_buf[0:msg_bits])
        {
            #pragma acc host_data use_device(rp, e1p, e2p)
            {
                rngongpu_fill_normal(rp,  (int)total,    1.0, noise_seed, 0xC0C0);
                rngongpu_fill_normal(e1p, (int)total,    1.0, noise_seed, 0xD0D0);
                rngongpu_fill_uniform_u32(e2p, (int)msg_bits, noise_seed, 0xE0E0);
            }

            #pragma acc parallel loop async(1) \
                present(rp[0:total], e1p[0:total], r_buf[0:total], e1_buf[0:total])
            for (size_t idx = 0; idx < total; ++idx) {
                double vr = rp[idx], ve = e1p[idx];
                if (vr >  6.0) vr =  6.0; if (vr < -6.0) vr = -6.0;
                if (ve >  6.0) ve =  6.0; if (ve < -6.0) ve = -6.0;
                r_buf[idx]  = (int32_t)(vr >= 0.0 ? vr + 0.5 : vr - 0.5);
                e1_buf[idx] = (int32_t)(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
            }

            #pragma acc parallel loop async(2) \
                present(e2p[0:msg_bits], e2_buf[0:msg_bits])
            for (uint32_t j = 0; j < msg_bits; ++j)
                e2_buf[j] = (int32_t)(e2p[j] % 3) - 1;

            #pragma acc wait(1,2)
        }
    }
}

/* ================================================================== */
/*  GenerateA_GPU_rngongpu — unchanged                                  */
/* ================================================================== */
void GenerateA_GPU_rngongpu(uint64_t rho_seed, uint32_t n, uint32_t q, int32_t* A_flat)
{
    size_t total = (size_t)n * n;
    #pragma acc parallel loop gang vector collapse(2) \
        firstprivate(rho_seed, n, q) present(A_flat[0:total])
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            A_flat[(size_t)i * n + j] = A_elem(rho_seed, (int)i, (int)j, (int)n, (int)q);
}

/* ================================================================== */
/*  KeyGen_GPU_rngongpu_Aflat                                           */
/*                                                                     */
/*  FIX 4: accepts pre-allocated persistent temp buffers.              */
/*  sraw_p (6n uint32_t), etmp_p (n double), ebuf_p (n int32_t)       */
/*  all device-resident from caller's persistent acc enter data.       */
/* ================================================================== */
void KeyGen_GPU_rngongpu_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                               int32_t* A_flat,
                               int32_t* s_out, int32_t* t_out,
                               uint32_t* sraw_p,  /* persistent: 6n uint32_t */
                               double*   etmp_p,  /* persistent: n double    */
                               int32_t*  ebuf_p)  /* persistent: n int32_t   */
{
    #pragma acc data present(A_flat[0:(size_t)n*n], s_out[0:n], t_out[0:n], \
                             sraw_p[0:6*n], etmp_p[0:n], ebuf_p[0:n])
    {
        /* Phase 1: secret s */
        #pragma acc host_data use_device(sraw_p)
        { rngongpu_fill_uniform_u32(sraw_p, (int)(6 * n), key_seed, 0xA0A0); }

        #pragma acc parallel loop present(sraw_p[0:6*n], s_out[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t base = 6 * i;
            int32_t a = (int32_t)(sraw_p[base+0]%2)+(int32_t)(sraw_p[base+1]%2)+(int32_t)(sraw_p[base+2]%2);
            int32_t b = (int32_t)(sraw_p[base+3]%2)+(int32_t)(sraw_p[base+4]%2)+(int32_t)(sraw_p[base+5]%2);
            s_out[i] = a - b;
        }

        /* Phase 2: Gaussian e */
        #pragma acc host_data use_device(etmp_p)
        { rngongpu_fill_normal(etmp_p, (int)n, 1, key_seed, 0xB0B0); }

        #pragma acc parallel loop present(etmp_p[0:n], ebuf_p[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            double ve = etmp_p[i];
            if (ve >  6.0) ve =  6.0;
            if (ve < -6.0) ve = -6.0;
            ebuf_p[i] = (int32_t)(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
        }

        /* Phase 3: t = A·s + e */
        #pragma acc parallel loop gang vector firstprivate(n, q) \
            present(A_flat[0:(size_t)n*n], s_out[0:n], ebuf_p[0:n], t_out[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t j = 0; j < n; ++j)
                acc += (long long)A_flat[(size_t)i * n + j] * (long long)s_out[j];
            long long v = (acc + (long long)ebuf_p[i]) % (long long)q;
            if (v < 0) v += q;
            t_out[i] = (int32_t)v;
        }
    }
}

/* ================================================================== */
/*  BatchEncrypt_GPU_Aflat — unchanged                                  */
/* ================================================================== */
void BatchEncrypt_GPU_Aflat(uint32_t n, uint32_t q,
                            int32_t* A_flat, int32_t* t_ptr,
                            int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                            int32_t* ptxt_ptr,
                            int32_t* c_out, uint32_t msg_bits)
{
    size_t total_r = (size_t)msg_bits * n;
    size_t total_c = (size_t)msg_bits * (n + 1);
    size_t total_A = (size_t)n * n;

    #pragma acc parallel loop gang firstprivate(n, q, msg_bits) \
        present(A_flat[0:total_A], t_ptr[0:n],                  \
                r_buf[0:total_r], e1_buf[0:total_r],            \
                e2_buf[0:msg_bits], ptxt_ptr[0:msg_bits],       \
                c_out[0:total_c])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t r_off = (size_t)j * n;
        size_t c_off = (size_t)j * (n + 1);

        #pragma acc loop vector
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t l = 0; l < n; ++l)
                acc += (long long)A_flat[(size_t)l * n + i] * (long long)r_buf[r_off + l];
            long long u = (acc + (long long)e1_buf[r_off + i]) % (long long)q;
            if (u < 0) u += q;
            c_out[c_off + i] = (int32_t)u;
        }

        long long dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += (long long)t_ptr[i] * (long long)r_buf[r_off + i];

        long long vv = (dot + (long long)e2_buf[j] + (long long)ptxt_ptr[j]) % (long long)q;
        if (vv < 0) vv += q;
        c_out[c_off + n] = (int32_t)vv;
    }
}

/* ================================================================== */
/*  KeyGen_GPU_rngongpu — on-the-fly A variant (local alloc, unchanged) */
/* ================================================================== */
void KeyGen_GPU_rngongpu(uint64_t key_seed, uint32_t n, uint32_t q,
                         int32_t* s_out, int32_t* t_out)
{
    std::vector<uint32_t> s_raw(6 * n);
    std::vector<double>   e_tmp(n);
    std::vector<int32_t>  e_buf(n);
    uint32_t* sraw_ptr = s_raw.data();
    double*   ep_tmp   = e_tmp.data();
    int32_t*  ebuf_ptr = e_buf.data();

    #pragma acc data create(sraw_ptr[0:6*n], ep_tmp[0:n], ebuf_ptr[0:n]) \
                    copyout(s_out[0:n], t_out[0:n])
    {
        #pragma acc host_data use_device(sraw_ptr)
        { rngongpu_fill_uniform_u32(sraw_ptr, (int)(6*n), key_seed, 0xA0A0); }

        #pragma acc parallel loop present(sraw_ptr[0:6*n], s_out[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t base = 6*i;
            int32_t a = (int32_t)(sraw_ptr[base+0]%2)+(int32_t)(sraw_ptr[base+1]%2)+(int32_t)(sraw_ptr[base+2]%2);
            int32_t b = (int32_t)(sraw_ptr[base+3]%2)+(int32_t)(sraw_ptr[base+4]%2)+(int32_t)(sraw_ptr[base+5]%2);
            s_out[i] = a - b;
        }

        #pragma acc host_data use_device(ep_tmp)
        { rngongpu_fill_normal(ep_tmp, (int)n, 1, key_seed, 0xB0B0); }

        #pragma acc parallel loop present(ep_tmp[0:n], ebuf_ptr[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            double ve = ep_tmp[i];
            if (ve > 6.0) ve = 6.0; if (ve < -6.0) ve = -6.0;
            ebuf_ptr[i] = (int32_t)(ve >= 0.0 ? ve + 0.5 : ve - 0.5);
        }

        #pragma acc parallel loop gang vector firstprivate(key_seed, n, q) \
            present(s_out[0:n], ebuf_ptr[0:n], t_out[0:n])
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t j = 0; j < n; ++j)
                acc += (long long)A_elem(key_seed, i, j, n, q) * (long long)s_out[j];
            long long v = (acc + (long long)ebuf_ptr[i]) % (long long)q;
            if (v < 0) v += q;
            t_out[i] = (int32_t)v;
        }
    }
}

/* ================================================================== */
/*  KeyGen_GPU — pure xorshift variant (unchanged)                     */
/* ================================================================== */
void KeyGen_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* s_out, int32_t* t_out)
{
    #pragma acc parallel loop firstprivate(key_seed, n) copyout(s_out[0:n])
    for (uint32_t i = 0; i < n; ++i) {
        uint64_t rng = key_seed ^ (0xA0A0A0A000000000ULL + i);
        s_out[i] = device_sample_cbd(&rng, 3);
    }

    #pragma acc parallel loop gang vector firstprivate(key_seed, n, q) \
        copyin(s_out[0:n]) copyout(t_out[0:n])
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += (long long)A_elem(key_seed, i, j, n, q) * (long long)s_out[j];
        uint64_t rng_e = key_seed ^ (0xB0B0B0B000000000ULL + i);
        int32_t ei = device_sample_gaussian(&rng_e, 2.3);
        long long v = (acc + ei) % (long long)q;
        if (v < 0) v += q;
        t_out[i] = (int32_t)v;
    }
}

/* ================================================================== */
/*  GPU_GenerateNoise — pure xorshift variant (unchanged)              */
/* ================================================================== */
void GPU_GenerateNoise(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                       int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf)
{
    size_t total = (size_t)msg_bits * n;

    #pragma acc parallel loop gang collapse(2) firstprivate(noise_seed, n, msg_bits) \
        copyout(r_buf[0:total], e1_buf[0:total])
    for (uint32_t j = 0; j < msg_bits; ++j)
        for (uint32_t i = 0; i < n; ++i) {
            uint64_t idx    = (uint64_t)j * n + i;
            uint64_t rng_r  = noise_seed ^ (0xC0C0C0C000000000ULL + idx);
            r_buf[idx]  = device_sample_gaussian(&rng_r, 1.0);
            uint64_t rng_e1 = noise_seed ^ (0xD0D0D0D000000000ULL + idx);
            e1_buf[idx] = device_sample_gaussian(&rng_e1, 1.0);
        }

    #pragma acc parallel loop firstprivate(noise_seed, msg_bits) copyout(e2_buf[0:msg_bits])
    for (uint32_t j = 0; j < msg_bits; ++j) {
        uint64_t rng_e2 = noise_seed ^ (0xE0E0E0E000000000ULL + j);
        e2_buf[j] = xorshift_uniform(&rng_e2, 2) - 1;
    }
}

/* ================================================================== */
/*  BatchEncrypt_GPU — on-the-fly A variant (unchanged)                */
/* ================================================================== */
void BatchEncrypt_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* t_ptr, int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                      int32_t* ptxt_ptr, int32_t* c_out, uint32_t msg_bits)
{
    size_t total_r = (size_t)msg_bits * n;
    size_t total_c = (size_t)msg_bits * (n + 1);

    #pragma acc parallel loop gang firstprivate(key_seed, n, q, msg_bits) \
        copyin(t_ptr[0:n], r_buf[0:total_r], e1_buf[0:total_r],  \
               e2_buf[0:msg_bits], ptxt_ptr[0:msg_bits])          \
        copyout(c_out[0:total_c])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t r_off = (size_t)j * n;
        size_t c_off = (size_t)j * (n + 1);

        #pragma acc loop vector
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t l = 0; l < n; ++l)
                acc += (long long)A_elem(key_seed, l, i, n, q) * (long long)r_buf[r_off + l];
            long long u = (acc + (long long)e1_buf[r_off + i]) % (long long)q;
            if (u < 0) u += q;
            c_out[c_off + i] = (int32_t)u;
        }

        long long dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += (long long)t_ptr[i] * (long long)r_buf[r_off + i];

        long long vv = (dot + (long long)e2_buf[j] + (long long)ptxt_ptr[j]) % (long long)q;
        if (vv < 0) vv += q;
        c_out[c_off + n] = (int32_t)vv;
    }
}

/* ================================================================== */
/*  BatchDecrypt_GPU — unchanged                                        */
/* ================================================================== */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits)
{
    size_t total_c = (size_t)msg_bits * (n + 1);

    #pragma acc parallel loop gang firstprivate(n, q, msg_bits) \
        present(s_ptr[0:n], c_in[0:total_c], decrypt_out[0:msg_bits])
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t c_off = (size_t)j * (n + 1);
        long long dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += (long long)s_ptr[i] * (long long)c_in[c_off + i];

        long long mu = ((long long)c_in[c_off + n] - dot) % (long long)q;
        if (mu < 0) mu += q;
        int32_t bound = (int32_t)q / 4;
        decrypt_out[j] = ((int32_t)mu <= bound || (int32_t)mu >= (int32_t)q - bound)
                           ? 0 : (int32_t)(q / 2);
    }
}

/* ================================================================== */
/*  CPU reference implementations (kept for debug Decrypt check)       */
/* ================================================================== */
void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>>& A,
            std::vector<int32_t>& s_k, std::vector<int32_t>& t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    std::vector<int32_t> Aflat = flatten_matrix(A);
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, 0);
    std::vector<int32_t> z(256);
    for (int i = 0; i < 256; ++i) z[i] = getRandomInt(0, 1);
    s_k = concat(s, z);
    int32_t* Af = Aflat.data(); int32_t* sp = s.data();
    int32_t* ep = e.data();     int32_t* pp = prod.data();
    #pragma acc parallel loop gang vector firstprivate(n, q) \
        copyin(Af[0:(size_t)n*n], sp[0:n], ep[0:n]) copyout(pp[0:n])
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += (long long)Af[i*n+j] * (long long)sp[j];
        long long v = (acc + (long long)ep[i]) % (long long)q;
        if (v < 0) v += q;
        pp[i] = (int32_t)v;
    }
    t = prod;
}

void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t>& t,
             std::vector<int32_t>& u, int32_t& v_i, uint32_t plaintext_i,
             std::vector<int32_t>& r, std::vector<int32_t>& e1, int32_t& e2,
             const std::vector<std::vector<int32_t>>& AT)
{
    std::vector<int32_t> ATflat = flatten_matrix(AT);
    u.resize(n);
    int32_t* ATf = ATflat.data(); int32_t* rp = r.data();
    int32_t* e1p = e1.data();     int32_t* up = u.data();
    int32_t* tp  = t.data();
    #pragma acc parallel loop gang vector firstprivate(n, q) \
        copyin(ATf[0:(size_t)n*n], rp[0:n], e1p[0:n]) copyout(up[0:n])
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += (long long)ATf[i*n+j] * (long long)rp[j];
        long long val = (acc + (long long)e1p[i]) % (long long)q;
        if (val < 0) val += q;
        up[i] = (int32_t)val;
    }
    long long dot = 0;
    #pragma acc parallel loop reduction(+:dot) firstprivate(n) copyin(tp[0:n], rp[0:n])
    for (uint32_t i = 0; i < n; ++i) dot += (long long)tp[i] * (long long)rp[i];
    long long res = dot % (long long)q; if (res < 0) res += q;
    long long vv = (res + (long long)e2 + (long long)plaintext_i) % (long long)q;
    if (vv < 0) vv += q;
    v_i = (int32_t)vv;
}

void Decrypt(int32_t v_i, const std::vector<int32_t>& u,
             const std::vector<int32_t>& s_k, uint32_t q, int32_t& decrypt_i)
{
    const size_t n = u.size();
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n);
    long long dot = 0;
    const int32_t* up = u.data(); const int32_t* sp = s.data();
    #pragma acc parallel loop reduction(+:dot) firstprivate(n) copyin(up[0:n], sp[0:n])
    for (size_t i = 0; i < n; ++i) dot += (long long)sp[i] * (long long)up[i];
    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0) r += q;
    int32_t mu = (int32_t)r;
    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)(q / 2);
}
