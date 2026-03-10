/*  pke.cpp  –  LWE public-key encryption (GPU via OpenMP target offloading)
 *
 *  GPU functions use xorshift64* to generate ALL cryptographic material
 *  on-device.  Matrix A is NEVER stored — each A[i][j] is recomputed
 *  from key_seed on the fly.
 *
 *  nvc++ 24.5 workarounds:
 *    1. Device functions embedded directly in this file (not from header)
 *       because nvc++ ICEs on #pragma omp declare target inline from headers
 *    2. No nested #pragma omp parallel for inside target teams distribute
 *       — split into separate target regions or use flat parallelism
 *    3. BatchEncrypt u-phase uses collapse(2), v-phase is separate target
 *
 *  Compile:  nvc++ -std=c++20 -O2 -mp=gpu -gpu=cc80,managed -Minfo=mp
 */
#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <cmath>
#include "pke.h"
#include "utils.h"

/* ================================================================== */
/*  Device functions — embedded in this TU for nvc++ compatibility     */
/*  (nvc++ 24.5 ICEs when these come from an included header)          */
/* ================================================================== */
#pragma omp begin declare target

static uint64_t dev_xorshift64star(uint64_t *state)
{
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

static int32_t dev_xorshift_uniform(uint64_t *state, int32_t max_val)
{
    return (int32_t)(dev_xorshift64star(state) % (uint64_t)(max_val + 1));
}

static int32_t dev_sample_cbd(uint64_t *state, int eta)
{
    int32_t a = 0, b = 0;
    for (int i = 0; i < eta; ++i) {
        a += dev_xorshift_uniform(state, 1);
        b += dev_xorshift_uniform(state, 1);
    }
    return a - b;
}

static int32_t dev_sample_gaussian(uint64_t *state, double sigma)
{
    double u1 = ((double)(dev_xorshift64star(state) % 1000000000ULL) + 1.0) / 1000000001.0;
    double u2 = ((double)(dev_xorshift64star(state) % 1000000000ULL) + 1.0) / 1000000001.0;

    double z = sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
    double val = sigma * z;

    double bnd = 6.0 * sigma;
    if (val >  bnd) val =  bnd;
    if (val < -bnd) val = -bnd;

    return (int32_t)(val >= 0.0 ? val + 0.5 : val - 0.5);
}

static int32_t dev_A_elem(uint64_t key_seed, int i, int j, int n, int q)
{
    uint64_t st = key_seed
                ^ ((uint64_t)i * 0x9E3779B97F4A7C15ULL)
                ^ ((uint64_t)j * 0x517CC1B727220A95ULL);
    return dev_xorshift_uniform(&st, q - 1);
}

#pragma omp end declare target

/* ================================================================== */
/*  KeyGen_GPU  –  s, t generated entirely on GPU                      */
/* ================================================================== */
void KeyGen_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* s_out, int32_t* t_out)
{
    /* Phase 1: generate secret vector s ~ CBD(eta=3)
       Simple flat parallel loop — no nesting issues */
    #pragma omp target teams distribute parallel for
    for (uint32_t i = 0; i < n; ++i) {
        uint64_t rng = key_seed ^ (0xA0A0A0A000000000ULL + i);
        s_out[i] = dev_sample_cbd(&rng, 3);
    }

    /* Phase 2: t = A · s + e  (mod q)
       Each thread computes one output t[i] with a sequential inner loop */
    #pragma omp target teams distribute parallel for
    for (uint32_t i = 0; i < n; ++i) {
        long long acc = 0;
        for (uint32_t j = 0; j < n; ++j)
            acc += (long long)dev_A_elem(key_seed, i, j, n, q) * (long long)s_out[j];

        uint64_t rng_e = key_seed ^ (0xB0B0B0B000000000ULL + i);
        int32_t ei = dev_sample_gaussian(&rng_e, 2.3);

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
    uint32_t total = msg_bits * n;

    #pragma omp target teams distribute parallel for
    for (uint32_t idx = 0; idx < total; ++idx) {
        uint64_t rng_r = noise_seed ^ (0xC0C0C0C000000000ULL + idx);
        r_buf[idx] = dev_sample_gaussian(&rng_r, 1.0);

        uint64_t rng_e1 = noise_seed ^ (0xD0D0D0D000000000ULL + idx);
        e1_buf[idx] = dev_sample_gaussian(&rng_e1, 1.0);
    }

    #pragma omp target teams distribute parallel for
    for (uint32_t j = 0; j < msg_bits; ++j) {
        uint64_t rng_e2 = noise_seed ^ (0xE0E0E0E000000000ULL + j);
        e2_buf[j] = dev_xorshift_uniform(&rng_e2, 6) - 3;
    }
}

/* ================================================================== */
/*  BatchEncrypt_GPU  –  256 encryptions, on-the-fly A                 */
/*                                                                     */
/*  Split into two separate target regions to avoid nested parallelism */
/*  that causes nvc++ ICE:                                             */
/*    Region 1: u = AT · r + e1  (collapse(2) over j×i)               */
/*    Region 2: v = t · r + e2 + ptxt  (one thread per j)             */
/* ================================================================== */
void BatchEncrypt_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* t_ptr,
                      int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                      int32_t* ptxt_ptr,
                      int32_t* c_out, uint32_t msg_bits)
{
    /* Region 1: u[j][i] = Σ_l A[l][i] * r[j*n+l] + e1[j*n+i]
       Flat 2D parallelism: each thread computes one (j, i) output */
    #pragma omp target teams distribute parallel for collapse(2)
    for (uint32_t j = 0; j < msg_bits; ++j) {
        for (uint32_t i = 0; i < n; ++i) {
            size_t r_off = (size_t)j * n;
            size_t c_off = (size_t)j * (n + 1);

            long long acc = 0;
            for (uint32_t l = 0; l < n; ++l)
                acc += (long long)dev_A_elem(key_seed, l, i, n, q)
                     * (long long)r_buf[r_off + l];

            long long u = (acc + (long long)e1_buf[r_off + i]) % (long long)q;
            if (u < 0) u += q;
            c_out[c_off + i] = (int32_t)u;
        }
    }

    /* Region 2: v[j] = t · r[j] + e2[j] + ptxt[j]  (mod q)
       Each thread computes one dot product sequentially */
    #pragma omp target teams distribute parallel for
    for (uint32_t j = 0; j < msg_bits; ++j) {
        size_t r_off = (size_t)j * n;
        size_t c_off = (size_t)j * (n + 1);

        long long dot = 0;
        for (uint32_t i = 0; i < n; ++i)
            dot += (long long)t_ptr[i] * (long long)r_buf[r_off + i];

        long long vv = (dot + (long long)e2_buf[j] + (long long)ptxt_ptr[j])
                       % (long long)q;
        if (vv < 0) vv += q;
        c_out[c_off + n] = (int32_t)vv;
    }
}

/* ================================================================== */
/*  BatchDecrypt_GPU  –  256 decryptions                               */
/*  Each thread handles one message bit with sequential reduction      */
/* ================================================================== */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits)
{
    #pragma omp target teams distribute parallel for
    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t c_off = (size_t)j * (n + 1);

        long long dot = 0;
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

    #pragma omp target teams distribute parallel for
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

    #pragma omp target teams distribute parallel for
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
    #pragma omp target teams distribute parallel for reduction(+:dot)
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

    #pragma omp target teams distribute parallel for reduction(+:dot)
    for (size_t i = 0; i < n; ++i)
        dot += (long long)sp[i] * (long long)up[i];

    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0) r += q;
    int32_t mu = (int32_t)r;
    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)(q / 2);
}
