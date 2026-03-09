/*  pke.cpp  –  LWE public-key encryption (GPU-accelerated)
 *
 *  GPU functions use xorshift64* to generate ALL cryptographic material
 *  on-device.  Matrix A is NEVER stored — each A[i][j] is recomputed
 *  from key_seed on the fly, eliminating ~640 MB of allocation and
 *  PCIe transfer for n=4096.
 *
 *  Compile:  nvc++ -acc -gpu=cc80,managed -Minfo=accel
 */
#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>
#include "pke.h"
#include "utils.h"
#include "xorshift.h"

/* ================================================================== */
/*  Device routine: deterministic A[i][j] from key_seed               */
/* ================================================================== */
#pragma acc routine seq
static inline int32_t A_elem(uint64_t key_seed, int i, int j, int n, int q)
{
    uint64_t st = key_seed
                ^ ((uint64_t)i * 0x9E3779B97F4A7C15ULL)
                ^ ((uint64_t)j * 0x517CC1B727220A95ULL);
    return xorshift_uniform(&st, q - 1);
}

/* ================================================================== */
/*  KeyGen_GPU  –  s, t generated entirely on GPU                      */
/* ================================================================== */
void KeyGen_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* s_out, int32_t* t_out)
{
    /* Phase 1: generate secret vector s ~ CBD(eta=3) */
    #pragma acc parallel loop firstprivate(key_seed, n)
    for (uint32_t i = 0; i < n; ++i) {
        uint64_t rng = key_seed ^ (0xA0A0A0A000000000ULL + i);
        s_out[i] = device_sample_cbd(&rng, 3);
    }

    /* Phase 2: t = A · s + e  (mod q)                                 */
    /*   A[i][j] regenerated on the fly — no storage needed             */
    #pragma acc parallel loop gang vector firstprivate(key_seed, n, q)
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
/*  GPU_GenerateNoise  –  r[msg_bits×n], e1[msg_bits×n], e2[msg_bits] */
/* ================================================================== */
void GPU_GenerateNoise(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                       int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf)
{
    #pragma acc parallel loop gang collapse(2) firstprivate(noise_seed, n, msg_bits)
    for (uint32_t j = 0; j < msg_bits; ++j) {
        for (uint32_t i = 0; i < n; ++i) {
            uint64_t idx = (uint64_t)j * n + i;

            uint64_t rng_r = noise_seed ^ (0xC0C0C0C000000000ULL + idx);
            r_buf[idx] = device_sample_gaussian(&rng_r, 1.0);

            uint64_t rng_e1 = noise_seed ^ (0xD0D0D0D000000000ULL + idx);
            e1_buf[idx] = device_sample_gaussian(&rng_e1, 1.0);
        }
    }

    #pragma acc parallel loop firstprivate(noise_seed, msg_bits)
    for (uint32_t j = 0; j < msg_bits; ++j) {
        uint64_t rng_e2 = noise_seed ^ (0xE0E0E0E000000000ULL + j);
        e2_buf[j] = xorshift_uniform(&rng_e2, 6) - 3;
    }
}

/* ================================================================== */
/*  BatchEncrypt_GPU  –  256 encryptions, 1 kernel, on-the-fly A      */
/* ================================================================== */
void BatchEncrypt_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* t_ptr,
                      int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                      int32_t* ptxt_ptr,
                      int32_t* c_out, uint32_t msg_bits)
{
    #pragma acc parallel loop gang firstprivate(key_seed, n, q, msg_bits)
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t r_off = (size_t)j * n;
        size_t c_off = (size_t)j * (n + 1);

        /* u = AT · r + e1  (AT[i][l] = A[l][i]) */
        #pragma acc loop vector
        for (uint32_t i = 0; i < n; ++i) {
            long long acc = 0;
            for (uint32_t l = 0; l < n; ++l)
                acc += (long long)A_elem(key_seed, l, i, n, q)
                     * (long long)r_buf[r_off + l];

            long long u = (acc + (long long)e1_buf[r_off + i]) % (long long)q;
            if (u < 0) u += q;
            c_out[c_off + i] = (int32_t)u;
        }

        /* v = t · r + e2 + ptxt  (mod q) */
        long long dot = 0;
        #pragma acc loop vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += (long long)t_ptr[i] * (long long)r_buf[r_off + i];

        long long vv = (dot + (long long)e2_buf[j] + (long long)ptxt_ptr[j])
                       % (long long)q;
        if (vv < 0) vv += q;
        c_out[c_off + n] = (int32_t)vv;
    }
}

/* ================================================================== */
/*  BatchDecrypt_GPU  –  256 decryptions, 1 kernel                     */
/* ================================================================== */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits)
{
    #pragma acc parallel loop gang firstprivate(n, q, msg_bits)
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
/*  Original CPU versions (kept for correctness reference)             */
/* ================================================================== */
void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A,
            std::vector<int32_t> &s_k, std::vector<int32_t> &t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    std::vector<int32_t> Aflat = flatten_matrix(A);
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, 0);

    std::vector<int32_t> z(256);
    for (int i = 0; i < 256; ++i) z[i] = getRandomInt(0, 1);
    s_k = concat(s, z);

    int32_t* Af = Aflat.data();
    int32_t* sp = s.data();
    int32_t* ep = e.data();
    int32_t* pp = prod.data();

    #pragma acc parallel loop gang vector firstprivate(n, q)
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += (long long)Af[i * n + j] * (long long)sp[j];
        long long v = (acc + (long long)ep[i]) % (long long)q;
        if (v < 0) v += q;
        pp[i] = (int32_t)v;
    }
    t = prod;
}

void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t,
             std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i,
             std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2,
             const std::vector<std::vector<int32_t>> &AT)
{
    std::vector<int32_t> ATflat = flatten_matrix(AT);
    u.resize(n);
    int32_t* ATf = ATflat.data();
    int32_t* rp  = r.data();
    int32_t* e1p = e1.data();
    int32_t* up  = u.data();

    #pragma acc parallel loop gang vector firstprivate(n, q)
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += (long long)ATf[i * n + j] * (long long)rp[j];
        long long val = (acc + (long long)e1p[i]) % (long long)q;
        if (val < 0) val += q;
        up[i] = (int32_t)val;
    }

    long long dot = 0;
    int32_t* tp = t.data();
    #pragma acc parallel loop reduction(+:dot) firstprivate(n)
    for (uint32_t i = 0; i < n; ++i)
        dot += (long long)tp[i] * (long long)rp[i];

    long long res = dot % (long long)q;
    if (res < 0) res += q;
    long long vv = (res + (long long)e2 + (long long)plaintext_i) % (long long)q;
    if (vv < 0) vv += q;
    v_i = (int32_t)vv;
}

void Decrypt(int32_t v_i, const std::vector<int32_t> &u,
             const std::vector<int32_t> &s_k, uint32_t q, int32_t &decrypt_i)
{
    const size_t n = u.size();
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n);
    long long dot = 0;
    const int32_t* up = u.data();
    const int32_t* sp = s.data();

    #pragma acc parallel loop reduction(+:dot) firstprivate(n)
    for (size_t i = 0; i < n; ++i)
        dot += (long long)sp[i] * (long long)up[i];

    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0) r += q;
    int32_t mu = (int32_t)r;
    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)(q / 2);
}
