#pragma once
#include <vector>
#include <cstdint>
#include "utils.h"

/* ---- Original CPU functions (kept for reference) ---- */
void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A,
            std::vector<int32_t> &s_k, std::vector<int32_t> &t);

void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u,
             int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r,
             std::vector<int32_t> &e1, int32_t &e2,
             const std::vector<std::vector<int32_t>> &AT);

void Decrypt(int32_t v_i, const std::vector<int32_t> &u,
             const std::vector<int32_t> &s_k, uint32_t q, int32_t &decrypt_i);

/* ---- GPU-accelerated functions (xorshift on-device) ---- */

/* KeyGen: generates s and t = A·s + e entirely on GPU.
   Matrix A is NEVER stored — computed on-the-fly from key_seed. */
void KeyGen_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* s_out, int32_t* t_out);

/* Generate all noise (r, e1, e2) on GPU from a seed */
void GPU_GenerateNoise(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                       int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf);

/* Batch encrypt: 256 encryptions in 1 kernel launch.
   A regenerated on-the-fly from key_seed (AT[i][l] = A[l][i]). */
void BatchEncrypt_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* t_ptr,
                      int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                      int32_t* ptxt_ptr,
                      int32_t* c_out, uint32_t msg_bits);

/* Batch decrypt: 256 decryptions in 1 kernel launch. */
void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits);