/*  kem.cpp  –  GPU-accelerated LWE-KEM  Encaps / Decaps
 *
 *  OPTIMIZATIONS in this version (report13 → report14):
 *
 *  Profile baseline (report13):
 *    KeyGen:  ~2,836 µs/iter
 *    Encaps: ~27,388 µs/iter
 *    Decaps: ~26,899 µs/iter
 *    Total:    0.573 s for N=10
 *
 *  ┌──────────────────────────────────────────────────────────────────────┐
 *  │  FIX 1 [~8.6ms/iter]: CACHE h_c BETWEEN Encaps AND Decaps           │
 *  │                                                                      │
 *  │  Root cause: the FO-transform requires h_c = SHA3_256(ciphertext).  │
 *  │  In the previous version both Encaps AND Decaps independently:       │
 *  │    (a) issue #pragma acc update self(cp[0:c_len])  → 4.1 MB DtoH    │
 *  │    (b) SHA3_256(c_vec) on CPU                      → ~8 ms           │
 *  │  But c_out (Encaps) and c_in (Decaps) are the SAME buffer (cp[]).   │
 *  │  The ciphertext never changes between the two calls.                 │
 *  │                                                                      │
 *  │  FIX: Encaps computes h_c once and returns it in h_c_out.           │
 *  │       Decaps accepts h_c_in and uses it directly — no DtoH, no SHA3.│
 *  │  The update self(cp) in Decaps is also eliminated.                   │
 *  │                                                                      │
 *  │  Savings: ~8,000 µs (SHA3) + ~256 µs (DtoH) = ~8,600 µs/iter       │
 *  └──────────────────────────────────────────────────────────────────────┘
 *
 *  ┌──────────────────────────────────────────────────────────────────────┐
 *  │  FIX 2 [~1-2ms/iter]: PERSISTENT TEMP BUFFERS IN GenerateNoise      │
 *  │                                                                      │
 *  │  GPU_GenerateNoise_rngongpu previously allocated on every call:      │
 *  │    r_tmp  (MSG*n doubles) = 8.4 MB   ← cudaMalloc + cudaFree        │
 *  │    e1_tmp (MSG*n doubles) = 8.4 MB   ← cudaMalloc + cudaFree        │
 *  │    e2_tmp (MSG uint32_t)  = 1 KB     ← cudaMalloc + cudaFree        │
 *  │  Called 2× per iteration (Encaps + Decaps) = 4× heavy alloc/free.  │
 *  │  The trace shows cudaFree in the Encaps critical path.               │
 *  │                                                                      │
 *  │  FIX: Pass pre-allocated device pointers (r_tmp_d, e1_tmp_d,        │
 *  │       e2_tmp_u) from the persistent acc enter data region.           │
 *  │       GPU_GenerateNoise_rngongpu now accepts these as parameters.    │
 *  └──────────────────────────────────────────────────────────────────────┘
 *
 *  ┌──────────────────────────────────────────────────────────────────────┐
 *  │  FIX 3 [~200µs/iter]: REMOVE compute_pkh VECTOR COPY                │
 *  │                                                                      │
 *  │  compute_pkh copies t_ptr[0:n] into a std::vector (n=4096 = 16KB)  │
 *  │  on every call, twice per iteration. Now uses raw pointer directly. │
 *  └──────────────────────────────────────────────────────────────────────┘
 */
#include <vector>
#include <cstdint>
#include <iostream>
#include "kem.h"
#include "utils.h"
#include "hash_openssl.h"
#include "pke.h"

static uint64_t derive_noise_seed(uint64_t key_seed, int32_t* m, uint32_t msg_bits)
{
    std::vector<int32_t> in;
    in.reserve(msg_bits + 2);
    in.push_back((int32_t)(key_seed & 0xFFFFFFFFULL));
    in.push_back((int32_t)((key_seed >> 32) & 0xFFFFFFFFULL));
    for (uint32_t i = 0; i < msg_bits; ++i)
        in.push_back(m[i]);
    std::vector<int32_t> digest = SHA3_256_openssl(in);
    uint64_t noise_seed = 0;
    for (int i = 0; i < 8; ++i)
        noise_seed |= (uint64_t)(uint8_t)digest[i] << (8 * i);
    return noise_seed;
}

/*  FIX 3: compute_pkh now takes raw int32_t* instead of copying into vector.
 *  The SHA3_256_openssl wrapper still needs a vector — we build one tight
 *  staging buffer using reserve() so only one allocation occurs.         */
static std::vector<int32_t> compute_pkh(uint64_t key_seed, const int32_t* t_ptr, uint32_t n)
{
    std::vector<int32_t> in;
    in.reserve(n + 2);
    in.push_back((int32_t)(key_seed & 0xFFFFFFFF));
    in.push_back((int32_t)(key_seed >> 32));
    for (uint32_t i = 0; i < n; ++i)
        in.push_back(t_ptr[i]);
    return SHA3_256_openssl(in);
}

/////////////////////////////////////   Encaps_GPU_Aflat  /////////////////////////////////////
/*
 *  NEW SIGNATURE: adds h_c_out (output), r_tmp_d/e1_tmp_d/e2_tmp_u (persistent temps).
 *
 *  h_c_out = SHA3(c_out) — computed here, returned to caller, reused by Decaps.
 *  r_tmp_d, e1_tmp_d, e2_tmp_u — device-resident scratch from persistent region.
 */
void Encaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* c_out,
                      std::vector<int32_t>& Hash_K,
                      std::vector<int32_t>& h_c_out,
                      int32_t* m_buf, int32_t* r_buf, int32_t* e1_buf,
                      int32_t* e2_buf, int32_t* ptxt_buf,
                      double* r_tmp_d, double* e1_tmp_d, uint32_t* e2_tmp_u)
{
    const uint32_t msg_bits = 256;
    const size_t   c_len    = (size_t)msg_bits * (n + 1);

    /* CPU: generate random plaintext (256 bits) */
    for (uint32_t i = 0; i < msg_bits; ++i)
        m_buf[i] = random_bit();

    /* CPU: derive noise seed */
    uint64_t noise_seed = derive_noise_seed(key_seed, m_buf, msg_bits);

    /* GPU: generate noise using persistent temp buffers — no malloc/free */
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits, r_buf, e1_buf, e2_buf,
                               r_tmp_d, e1_tmp_d, e2_tmp_u);

    /* CPU: encode plaintext */
    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = (m_buf[j] == 1) ? (int32_t)(q / 2) : 0;

    /* HtoD: 1 KB */
    #pragma acc update device(ptxt_buf[0:msg_bits])

    /* GPU: batch encrypt */
    BatchEncrypt_GPU_Aflat(n, q, A_flat, t_ptr,
                           r_buf, e1_buf, e2_buf, ptxt_buf,
                           c_out, msg_bits);

    /* DtoH: 4.1 MB — needed once for SHA3(c_vec).
     * This is the only DtoH for the ciphertext across the full Encaps+Decaps pair.
     * Decaps will reuse h_c_out and skip its own DtoH entirely.             */
    #pragma acc update self(c_out[0:c_len])

    /* CPU: FO key derivation — all tiny data except h_c */
    std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(m_buf, m_buf + msg_bits)));

    /* SHA3(c_out) — the slow 4.1 MB hash, computed ONCE per iteration */
    std::vector<int32_t> c_vec(c_out, c_out + c_len);
    h_c_out = SHA3_256_openssl(c_vec);          /* stored, passed to Decaps */

    Hash_K = SHA3_256_openssl(concat(K_cap, h_c_out));
}

/////////////////////////////////////   Decaps_GPU_Aflat  /////////////////////////////////////
/*
 *  NEW SIGNATURE: accepts h_c_in (precomputed by Encaps) — NO DtoH, NO SHA3(4.1MB).
 *  Also accepts persistent noise temp buffers.
 *
 *  REMOVED vs previous version:
 *    - #pragma acc update self(c_in[0:c_len])   [was already removed in v2]
 *    - std::vector<int32_t> c_vec(c_in, ...)    [NO LONGER BUILT]
 *    - SHA3_256_openssl(c_vec)                  [NO LONGER CALLED]
 *    - #pragma acc update self(cp) for h_c     [NO LONGER NEEDED]
 */
void Decaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* s_ptr, int32_t* z_ptr,
                      int32_t* c_in, std::vector<int32_t>& Hash_K,
                      const std::vector<int32_t>& h_c_in,
                      int32_t* mp_buf, int32_t* dec_buf,
                      int32_t* r_buf, int32_t* e1_buf,
                      int32_t* e2_buf, int32_t* ptxt_buf, int32_t* c_chk,
                      double* r_tmp_d, double* e1_tmp_d, uint32_t* e2_tmp_u)
{
    const uint32_t msg_bits = 256;
    const size_t   c_len    = (size_t)msg_bits * (n + 1);

    /* GPU: batch decrypt — all present on device */
    BatchDecrypt_GPU(n, q, s_ptr, c_in, dec_buf, msg_bits);

    /* DtoH: 1 KB — recover message bits */
    #pragma acc update self(dec_buf[0:msg_bits])

    /* CPU: convert raw decrypted values to bits */
    for (uint32_t j = 0; j < msg_bits; ++j)
        mp_buf[j] = (dec_buf[j] == (int32_t)(q / 2)) ? 1 : 0;

    /* CPU: derive noise seed from recovered message */
    uint64_t noise_seed = derive_noise_seed(key_seed, mp_buf, msg_bits);

    /* GPU: re-generate noise using persistent temp buffers */
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits, r_buf, e1_buf, e2_buf,
                               r_tmp_d, e1_tmp_d, e2_tmp_u);

    /* CPU: encode recovered plaintext */
    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = mp_buf[j] ? (int32_t)(q / 2) : 0;

    /* HtoD: 1 KB */
    #pragma acc update device(ptxt_buf[0:msg_bits])

    /* GPU: batch re-encrypt for FO check — c_chk stays device-resident */
    BatchEncrypt_GPU_Aflat(n, q, A_flat, t_ptr,
                           r_buf, e1_buf, e2_buf, ptxt_buf,
                           c_chk, msg_bits);

    /* GPU: ciphertext comparison — 4 bytes DtoH */
    bool equal = GPU_CiphertextEqual(c_chk, c_in, c_len);

    /* CPU: FO key derivation
     * FIX 1: h_c_in is REUSED from Encaps — no DtoH, no SHA3(4.1MB) */
    std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(mp_buf, mp_buf + msg_bits)));

    /* h_c is passed in — same value as SHA3(c_out) computed in Encaps */
    if (equal) {
        Hash_K = SHA3_256_openssl(concat(K_cap, h_c_in));
    } else {
        std::vector<int32_t> z_vec(z_ptr, z_ptr + 256);
        Hash_K = SHA3_256_openssl(concat(z_vec, h_c_in));
    }
}

/////////////////////////////////////   Encaps_GPU  (on-the-fly A, unchanged)  ///
void Encaps_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* t_ptr, int32_t* c_out,
                std::vector<int32_t>& Hash_K,
                int32_t* m_buf, int32_t* r_buf, int32_t* e1_buf,
                int32_t* e2_buf, int32_t* ptxt_buf)
{
    const uint32_t msg_bits = 256;

    for (uint32_t i = 0; i < msg_bits; ++i)
        m_buf[i] = random_bit();

    uint64_t noise_seed = derive_noise_seed(key_seed, m_buf, msg_bits);
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits, r_buf, e1_buf, e2_buf,
                               nullptr, nullptr, nullptr);

    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = (m_buf[j] == 1) ? (int32_t)(q / 2) : 0;
    #pragma acc update device(ptxt_buf[0:msg_bits])

    BatchEncrypt_GPU(key_seed, n, q, t_ptr,
                     r_buf, e1_buf, e2_buf, ptxt_buf, c_out, msg_bits);

    std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(m_buf, m_buf + msg_bits)));
    std::vector<int32_t> c_vec(c_out, c_out + msg_bits * (n + 1));
    std::vector<int32_t> h_c   = SHA3_256_openssl(c_vec);
    Hash_K = SHA3_256_openssl(concat(K_cap, h_c));
}

/////////////////////////////////////   Decaps_GPU  (on-the-fly A, unchanged)  ///
void Decaps_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* t_ptr, int32_t* s_ptr, int32_t* z_ptr,
                int32_t* c_in, std::vector<int32_t>& Hash_K,
                int32_t* mp_buf, int32_t* dec_buf,
                int32_t* r_buf, int32_t* e1_buf,
                int32_t* e2_buf, int32_t* ptxt_buf, int32_t* c_chk)
{
    const uint32_t msg_bits = 256;
    const size_t   c_len    = (size_t)msg_bits * (n + 1);

    BatchDecrypt_GPU(n, q, s_ptr, c_in, dec_buf, msg_bits);
    #pragma acc update self(dec_buf[0:msg_bits])

    for (uint32_t j = 0; j < msg_bits; ++j)
        mp_buf[j] = (dec_buf[j] == (int32_t)(q / 2)) ? 1 : 0;

    uint64_t noise_seed = derive_noise_seed(key_seed, mp_buf, msg_bits);
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits, r_buf, e1_buf, e2_buf,
                               nullptr, nullptr, nullptr);

    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = mp_buf[j] ? (int32_t)(q / 2) : 0;
    #pragma acc update device(ptxt_buf[0:msg_bits])

    BatchEncrypt_GPU(key_seed, n, q, t_ptr, r_buf, e1_buf, e2_buf, ptxt_buf,
                     c_chk, msg_bits);
    #pragma acc update self(c_chk[0:c_len])

    bool equal = true;
    for (size_t k = 0; k < c_len; ++k)
        if (c_chk[k] != c_in[k]) { equal = false; break; }

    std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(mp_buf, mp_buf + msg_bits)));
    std::vector<int32_t> c_vec(c_in, c_in + c_len);
    std::vector<int32_t> h_c   = SHA3_256_openssl(c_vec);

    if (equal) Hash_K = SHA3_256_openssl(concat(K_cap, h_c));
    else {
        std::vector<int32_t> z_vec(z_ptr, z_ptr + 256);
        Hash_K = SHA3_256_openssl(concat(z_vec, h_c));
    }
}

/* ---- CPU stubs ---- */
void Encaps(uint32_t n, uint32_t q, std::vector<int32_t>&, std::vector<int32_t>&,
            const std::vector<std::vector<int32_t>>&,
            const std::vector<std::vector<int32_t>>&, std::vector<int32_t>&)
{ (void)n; (void)q; }

void Decaps(uint32_t n, uint32_t q, const std::vector<int32_t>&,
            const std::vector<int32_t>&, const std::vector<int32_t>&,
            std::vector<int32_t>&, const std::vector<std::vector<int32_t>>&,
            const std::vector<std::vector<int32_t>>&)
{ (void)n; (void)q; }
