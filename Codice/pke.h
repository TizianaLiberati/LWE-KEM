/*  pke.h  –  GPU-accelerated LWE PKE/KEM primitives                  */
#pragma once
#include <vector>
#include <cstdint>

/* ---- Key generation ---- */
void KeyGen_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* s_out, int32_t* t_out);

void KeyGen_GPU_rngongpu(uint64_t key_seed, uint32_t n, uint32_t q,
                         int32_t* s_out, int32_t* t_out);

/* FIX 4: accepts persistent temp buffers — no malloc on hot path */
void KeyGen_GPU_rngongpu_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                               int32_t* A_flat,
                               int32_t* s_out, int32_t* t_out,
                               uint32_t* sraw_p,   /* persistent: 6n uint32_t */
                               double*   etmp_p,   /* persistent: n double    */
                               int32_t*  ebuf_p);  /* persistent: n int32_t   */

void KeyGen(uint32_t n, uint32_t q,
            std::vector<std::vector<int32_t>>& A,
            std::vector<int32_t>& s_k,
            std::vector<int32_t>& t);

/* ---- Matrix generation ---- */
void GenerateA_GPU_rngongpu(uint64_t rho_seed, uint32_t n, uint32_t q, int32_t* A_flat);

/* ---- Noise generation ---- */
void GPU_GenerateNoise(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                       int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf);

/* FIX 2: accepts nullable persistent device pointers for double/u32 intermediates */
void GPU_GenerateNoise_rngongpu(uint64_t noise_seed, uint32_t n, uint32_t msg_bits,
                                int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                                double*   r_tmp_d,   /* nullable */
                                double*   e1_tmp_d,  /* nullable */
                                uint32_t* e2_tmp_u); /* nullable */

/* ---- Batch encrypt/decrypt ---- */
void BatchEncrypt_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* t_ptr,
                      int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                      int32_t* ptxt_ptr, int32_t* c_out, uint32_t msg_bits);

void BatchEncrypt_GPU_Aflat(uint32_t n, uint32_t q,
                            int32_t* A_flat, int32_t* t_ptr,
                            int32_t* r_buf, int32_t* e1_buf, int32_t* e2_buf,
                            int32_t* ptxt_ptr, int32_t* c_out, uint32_t msg_bits);

void BatchDecrypt_GPU(uint32_t n, uint32_t q,
                      int32_t* s_ptr, int32_t* c_in,
                      int32_t* decrypt_out, uint32_t msg_bits);

/* ---- FO-transform helper ---- */
bool GPU_CiphertextEqual(const int32_t* c1, const int32_t* c2, size_t len);

/* ---- CPU reference implementations ---- */
void Encrypt(uint32_t n, uint32_t q,
             std::vector<int32_t>& t, std::vector<int32_t>& u,
             int32_t& v_i, uint32_t plaintext_i,
             std::vector<int32_t>& r, std::vector<int32_t>& e1, int32_t& e2,
             const std::vector<std::vector<int32_t>>& AT);

void Decrypt(int32_t v_i, const std::vector<int32_t>& u,
             const std::vector<int32_t>& s_k, uint32_t q, int32_t& decrypt_i);
