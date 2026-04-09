/*  kem.cpp  –  GPU-accelerated LWE-KEM  Encaps / Decaps  (v4)
 *
 *  v4 changes vs v3:
 *    • A_flat parameter type: int32_t* → uint16_t* (OPT-5).
 *    • r_buf / e1_buf parameters renamed r_buf_T / e1_buf_T and now hold
 *      the l-major (transposed) layout produced by GPU_GenerateNoise_rngongpu
 *      (OPT-8).
 *    • BatchEncrypt_GPU_Aflat signature updated to accept pool& for AT access
 *      (OPT-6).
 *    • No logic changes — correctness is preserved.
 */
#include <vector>
#include <cstdint>
#include <iostream>
#include "kem.h"
#include "utils.h"
#include "hash_openssl.h"
#include "pke.h"

/* ── Internal helpers ─────────────────────────────────────────────────────── */

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

static std::vector<int32_t> compute_pkh(uint64_t key_seed,
                                        int32_t* t_ptr, uint32_t n)
{
    std::vector<int32_t> in;
    in.reserve(n + 2);
    in.push_back((int32_t)(key_seed & 0xFFFFFFFF));
    in.push_back((int32_t)(key_seed >> 32));
    for (uint32_t i = 0; i < n; ++i)
        in.push_back(t_ptr[i]);
    return SHA3_256_openssl(in);
}

/* ══════════════════════════════════════════════════════════════════════════
 *  Encaps_GPU_Aflat
 *
 *  OPT-5: A_flat is uint16_t*.
 *  OPT-6/8: r_buf_T / e1_buf_T carry the l-major transposed noise;
 *            BatchEncrypt_GPU_Aflat reads pool.AT internally.
 * ══════════════════════════════════════════════════════════════════════════ */
void Encaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      uint16_t* A_flat, int32_t* t_ptr, int32_t* c_out,
                      std::vector<int32_t>& Hash_K,
                      int32_t* m_buf,
                      int32_t* r_buf_T, int32_t* e1_buf_T,
                      int32_t* e2_buf,  int32_t* ptxt_buf,
                      DevicePool& pool)
{
    const uint32_t msg_bits = 256;

    for (uint32_t i = 0; i < msg_bits; ++i)
        m_buf[i] = random_bit();

    uint64_t noise_seed = derive_noise_seed(key_seed, m_buf, msg_bits);

    /* GPU: fill r_buf_T / e1_buf_T in l-major layout */
    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits,
                               r_buf_T, e1_buf_T, e2_buf, pool);

    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = (m_buf[j] == 1) ? (int32_t)(q / 2) : 0;

    /* GPU: encrypt using pool.AT + transposed r (OPT-6 + OPT-8) */
    BatchEncrypt_GPU_Aflat(n, q, A_flat, t_ptr,
                           r_buf_T, e1_buf_T, e2_buf, ptxt_buf,
                           c_out, msg_bits, pool);

    std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(m_buf, m_buf + msg_bits)));
    std::vector<int32_t> c_vec(c_out, c_out + (size_t)msg_bits * (n + 1));
    std::vector<int32_t> h_c   = SHA3_256_openssl(c_vec);
    Hash_K = SHA3_256_openssl(concat(K_cap, h_c));
}

/* ══════════════════════════════════════════════════════════════════════════
 *  Decaps_GPU_Aflat
 *
 *  OPT-5: A_flat is uint16_t*.
 *  OPT-6/8: same transposed noise layout as Encaps.
 * ══════════════════════════════════════════════════════════════════════════ */
void Decaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      uint16_t* A_flat, int32_t* t_ptr,
                      int32_t* s_ptr,   int32_t* z_ptr,
                      int32_t* c_in,    std::vector<int32_t>& Hash_K,
                      int32_t* mp_buf,  int32_t* dec_buf,
                      int32_t* r_buf_T, int32_t* e1_buf_T,
                      int32_t* e2_buf,  int32_t* ptxt_buf,
                      int32_t* c_chk,
                      DevicePool& pool)
{
    const uint32_t msg_bits = 256;
    const size_t   c_len    = (size_t)msg_bits * (n + 1);

    BatchDecrypt_GPU(n, q, s_ptr, c_in, dec_buf, msg_bits);

    for (uint32_t j = 0; j < msg_bits; ++j)
        mp_buf[j] = (dec_buf[j] == (int32_t)(q / 2)) ? 1 : 0;

    uint64_t noise_seed = derive_noise_seed(key_seed, mp_buf, msg_bits);

    GPU_GenerateNoise_rngongpu(noise_seed, n, msg_bits,
                               r_buf_T, e1_buf_T, e2_buf, pool);

    for (uint32_t j = 0; j < msg_bits; ++j)
        ptxt_buf[j] = mp_buf[j] ? (int32_t)(q / 2) : 0;

    BatchEncrypt_GPU_Aflat(n, q, A_flat, t_ptr,
                           r_buf_T, e1_buf_T, e2_buf, ptxt_buf,
                           c_chk, msg_bits, pool);

    bool equal = true;
    for (size_t k = 0; k < c_len; ++k)
        if (c_chk[k] != c_in[k]) { equal = false; break; }

    std::vector<int32_t> pkh   = compute_pkh(key_seed, t_ptr, n);
    std::vector<int32_t> K_cap = SHA3_256_openssl(
        concat(pkh, std::vector<int32_t>(mp_buf, mp_buf + msg_bits)));
    std::vector<int32_t> c_vec(c_in, c_in + c_len);
    std::vector<int32_t> h_c   = SHA3_256_openssl(c_vec);

    if (equal) {
        Hash_K = SHA3_256_openssl(concat(K_cap, h_c));
    } else {
        std::vector<int32_t> z_vec(z_ptr, z_ptr + 256);
        Hash_K = SHA3_256_openssl(concat(z_vec, h_c));
    }
}

/* ── Stub CPU versions (no logic change) ─────────────────────────────────── */
void Encaps(uint32_t, uint32_t, std::vector<int32_t>&, std::vector<int32_t>&,
            const std::vector<std::vector<int32_t>>&,
            const std::vector<std::vector<int32_t>>&,
            std::vector<int32_t>&) {}

void Decaps(uint32_t, uint32_t,
            const std::vector<int32_t>&, const std::vector<int32_t>&,
            const std::vector<int32_t>&, std::vector<int32_t>&,
            const std::vector<std::vector<int32_t>>&,
            const std::vector<std::vector<int32_t>>&) {}
