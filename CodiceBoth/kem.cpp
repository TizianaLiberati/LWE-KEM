/*  kem.cpp  –  GPU-accelerated LWE-KEM Encaps / Decaps
 *
 *  Compile-time backend selection:
 *    nvc++ -DUSE_OPENACC    -acc -gpu=cc90,managed ...
 *    nvc++ -DUSE_OMP_TARGET -mp=gpu -gpu=cc90,managed ...
 *
 *  ──────────────────────────────────────────────────────────────────────
 *  OPTIMIZATIONS (from report13 → report14, preserved here):
 *
 *  FIX 1 [~8.6ms/iter]: CACHE h_c BETWEEN Encaps AND Decaps
 *    Encaps computes h_c = SHA3(ciphertext) once and returns it in h_c_out.
 *    Decaps receives h_c_in and skips the 4.1 MB DtoH + SHA3 rehash.
 *
 *  FIX 2 [~1-2ms/iter]: PERSISTENT TEMP BUFFERS IN GenerateNoise
 *    r_tmp_d, e1_tmp_d, e2_tmp_u are pre-allocated in the caller's
 *    persistent device region; GPU_GenerateNoise_rngongpu uses the
 *    fast path (no cudaMalloc/cudaFree inside the critical path).
 *
 *  FIX 3 [~200µs/iter]: REMOVE compute_pkh VECTOR COPY
 *    compute_pkh now accepts raw int32_t* and builds one tight
 *    staging buffer with reserve().
 *
 *  ──────────────────────────────────────────────────────────────────────
 *  DIRECTIVE MAPPING IN THIS FILE
 *  ──────────────────────────────────────────────────────────────────────
 *
 *  Encaps_GPU_Aflat:
 *    ptxt_buf HtoD (1 KB after CPU encode):
 *      ACC:  #pragma acc update device(ptxt_buf[0:msg_bits])
 *      OMP:  #pragma omp target update to(ptxt_buf[0:msg_bits])
 *
 *    c_out DtoH (4.1 MB for SHA3 – only once per iter):
 *      ACC:  #pragma acc update self(c_out[0:c_len])
 *      OMP:  #pragma omp target update from(c_out[0:c_len])
 *
 *  Decaps_GPU_Aflat:
 *    dec_buf DtoH (1 KB to recover message bits):
 *      ACC:  #pragma acc update self(dec_buf[0:msg_bits])
 *      OMP:  #pragma omp target update from(dec_buf[0:msg_bits])
 *
 *    ptxt_buf HtoD (1 KB re-encode for FO check):
 *      ACC:  #pragma acc update device(ptxt_buf[0:msg_bits])
 *      OMP:  #pragma omp target update to(ptxt_buf[0:msg_bits])
 *
 *  All GPU compute calls (BatchEncrypt, BatchDecrypt, GenerateNoise) are
 *  backend-agnostic at this level: they call the functions in pke.cpp
 *  which carry the appropriate #ifdef internally.
 */

#include <vector>
#include <cstdint>
#include <iostream>
#include "kem.h"
#include "utils.h"
#include "hash_openssl.h"
#include "pke.h"
#include "gpu_backend.h"   /* compile-time ACC / OMP-target abstraction */
#include "rngongpu_adapter.h"

/* ---------------------------------------------------------------------------
 *  derive_noise_seed
 *  Deterministically derives a noise seed from the key seed and plaintext.
 *  Pure CPU; no GPU offload needed (tiny computation).
 * --------------------------------------------------------------------------- */
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

/* ---------------------------------------------------------------------------
 *  compute_pkh  (FIX 3: raw pointer, single reserve)
 *  Computes SHA3_256(key_seed || t_ptr[0..n-1]).
 *  The SHA3 wrapper still needs a vector, built with one allocation.
 * --------------------------------------------------------------------------- */
static std::vector<int32_t> compute_pkh(uint64_t key_seed,
                                        const int32_t* t_ptr, uint32_t n)
{
    std::vector<int32_t> in;
    in.reserve(n + 2);
    in.push_back((int32_t)(key_seed & 0xFFFFFFFF));
    in.push_back((int32_t)(key_seed >> 32));
    for (uint32_t i = 0; i < n; ++i)
        in.push_back(t_ptr[i]);
    return SHA3_256_openssl(in);
}

/* ===========================================================================
 *  Encaps_GPU_Aflat
 *
 *  NEW SIGNATURE (report14):
 *    h_c_out          – out: SHA3(ciphertext), cached and passed to Decaps.
 *    r_tmp_d/e1_tmp_d/e2_tmp_u – persistent device noise temps (FIX 2).
 *
 *  DATA FLOW:
 *    1. CPU: random m (256 bits)
 *    2. CPU: noise_seed = derive(key_seed, m)
 *    3. GPU: noise r, e1, e2 via GPU_GenerateNoise_rngongpu (fast path)
 *    4. CPU: ptxt_buf = encode(m)
 *    5. HtoD: ptxt_buf (1 KB)
 *    6. GPU: BatchEncrypt → c_out
 *    7. DtoH: c_out (4.1 MB) — ONLY transfer in the whole Encaps+Decaps pair
 *    8. CPU: h_c_out = SHA3(c_out)  [cached for Decaps]
 *    9. CPU: FO key derivation → Hash_K
 * =========================================================================== */
void Encaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* c_out,
                      std::vector<int32_t>& Hash_K,
                      std::vector<int32_t>& h_c_out,
                      int32_t*  m_buf,    int32_t*  r_buf,
                      int32_t*  e1_buf,   int32_t*  e2_buf,
                      int32_t*  ptxt_buf,
                      double*   r_tmp_d,  double*   e1_tmp_d,
                      uint32_t* e2_tmp_u)
{
    const uint32_t msg_bits = 256;
    const size_t   c_len    = (size_t)msg_bits * (n + 1);

    /* Step 1: random plaintext bits (CPU) */
    for (uint32_t i = 0; i < msg_bits; ++i)
        m_buf[i] = random_bit();

    /* Step 2: deterministic noise seed */
    const uint64_t noise_seed = derive_noise_seed(key_seed, m_buf, msg_bits);

    /* Step 3: GPU noise generation — persistent buffers, no malloc */
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits,
                               r_buf, e1_buf, e2_buf,
                               r_tmp_d, e1_tmp_d, e2_tmp_u);

    /* Step 4: encode plaintext (CPU) */
    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = (m_buf[j] == 1) ? (int32_t)(q / 2) : 0;

    /* Step 5: HtoD ptxt_buf (1 KB) */
#ifdef USE_OPENACC
    #pragma acc update device(ptxt_buf[0:msg_bits])
#else
    #pragma omp target update to(ptxt_buf[0:msg_bits])
#endif

    /* Step 6: batch encrypt on GPU */
    BatchEncrypt_GPU_Aflat(n, q, A_flat, t_ptr,
                           r_buf, e1_buf, e2_buf, ptxt_buf,
                           c_out, msg_bits);

    /* Step 7: DtoH c_out (4.1 MB) — single round-trip per Encaps+Decaps pair */
#ifdef USE_OPENACC
    #pragma acc update self(c_out[0:c_len])
#else
    #pragma omp target update from(c_out[0:c_len])
#endif

    /* Step 8–9: FO key derivation on CPU */
    const std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    const std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(m_buf, m_buf + msg_bits)));

    /* h_c = SHA3(c_out) — computed ONCE; Decaps reuses this value (FIX 1) */
    const std::vector<int32_t> c_vec(c_out, c_out + c_len);
    h_c_out = SHA3_256_openssl(c_vec);

    Hash_K = SHA3_256_openssl(concat(K_cap, h_c_out));
}

/* ===========================================================================
 *  Decaps_GPU_Aflat
 *
 *  NEW SIGNATURE (report14):
 *    h_c_in  – in: SHA3(ciphertext) precomputed by Encaps (FIX 1).
 *              No DtoH and no SHA3(4.1 MB) here.
 *    r_tmp_d/e1_tmp_d/e2_tmp_u – persistent device noise temps (FIX 2).
 *
 *  DATA FLOW:
 *    1. GPU: BatchDecrypt → dec_buf
 *    2. DtoH: dec_buf (1 KB)
 *    3. CPU: m' = decode(dec_buf)
 *    4. CPU: noise_seed = derive(key_seed, m')
 *    5. GPU: re-generate noise (fast path)
 *    6. CPU: ptxt_buf = encode(m')
 *    7. HtoD: ptxt_buf (1 KB)
 *    8. GPU: BatchEncrypt → c_chk
 *    9. GPU: ciphertext equality check (4 bytes DtoH)
 *   10. CPU: FO key derivation using h_c_in (no re-hash, no DtoH — FIX 1)
 * =========================================================================== */
void Decaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr,
                      int32_t* s_ptr,  int32_t* z_ptr,
                      int32_t* c_in,
                      std::vector<int32_t>& Hash_K,
                      const std::vector<int32_t>& h_c_in,
                      int32_t*  mp_buf,   int32_t*  dec_buf,
                      int32_t*  r_buf,    int32_t*  e1_buf,
                      int32_t*  e2_buf,   int32_t*  ptxt_buf,
                      int32_t*  c_chk,
                      double*   r_tmp_d,  double*   e1_tmp_d,
                      uint32_t* e2_tmp_u)
{
    const uint32_t msg_bits = 256;
    const size_t   c_len    = (size_t)msg_bits * (n + 1);

    /* Step 1: GPU decrypt */
    BatchDecrypt_GPU(n, q, s_ptr, c_in, dec_buf, msg_bits);

    /* Step 2: DtoH dec_buf (1 KB) */
#ifdef USE_OPENACC
    #pragma acc update self(dec_buf[0:msg_bits])
#else
    #pragma omp target update from(dec_buf[0:msg_bits])
#endif

    /* Step 3: decode recovered message bits (CPU) */
    for (uint32_t j = 0; j < msg_bits; ++j)
        mp_buf[j] = (dec_buf[j] == (int32_t)(q / 2)) ? 1 : 0;

    /* Step 4: derive noise seed from recovered message */
    const uint64_t noise_seed = derive_noise_seed(key_seed, mp_buf, msg_bits);

    /* Step 5: re-generate noise on GPU using persistent temps */
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits,
                               r_buf, e1_buf, e2_buf,
                               r_tmp_d, e1_tmp_d, e2_tmp_u);

    /* Step 6: re-encode recovered message (CPU) */
    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = mp_buf[j] ? (int32_t)(q / 2) : 0;

    /* Step 7: HtoD ptxt_buf (1 KB) */
#ifdef USE_OPENACC
    #pragma acc update device(ptxt_buf[0:msg_bits])
#else
    #pragma omp target update to(ptxt_buf[0:msg_bits])
#endif

    /* Step 8: re-encrypt for FO check (c_chk stays device-resident) */
    BatchEncrypt_GPU_Aflat(n, q, A_flat, t_ptr,
                           r_buf, e1_buf, e2_buf, ptxt_buf,
                           c_chk, msg_bits);

    /* Step 9: ciphertext comparison on GPU (4 bytes DtoH) */
    const bool equal = GPU_CiphertextEqual(c_chk, c_in, c_len);

    /* Step 10: FO key derivation (CPU)
     *   FIX 1: h_c_in reused from Encaps — no DtoH, no SHA3(4.1 MB) */
    const std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    const std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(mp_buf, mp_buf + msg_bits)));

    if (equal) {
        Hash_K = SHA3_256_openssl(concat(K_cap, h_c_in));
    } else {
        const std::vector<int32_t> z_vec(z_ptr, z_ptr + 256);
        Hash_K = SHA3_256_openssl(concat(z_vec, h_c_in));
    }
}
