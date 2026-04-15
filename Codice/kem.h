/*  kem.h  –  KEM interface with h_c caching to eliminate duplicate SHA3(4.1MB)  */
#pragma once
#include <vector>
#include <cstdint>

/*
 *  OPTIMIZATION: Encaps now returns h_c = SHA3(c_out) alongside Hash_K.
 *  The caller passes h_c directly to Decaps, avoiding:
 *    - One acc update self(cp) → 4.1 MB DtoH
 *    - One SHA3(c_in = 4.1 MB) CPU call (~8 ms)
 *  Total savings: ~8.6 ms per iteration.
 *
 *  h_c is a 32-element int32_t vector (SHA3-256 output = 32 bytes).
 */

/* ---- GPU-accelerated variants with precomputed A_flat ---- */
void Encaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* c_out,
                      std::vector<int32_t>& Hash_K,
                      std::vector<int32_t>& h_c_out,   /* NEW: cached h_c */
                      int32_t* m_buf, int32_t* r_buf, int32_t* e1_buf,
                      int32_t* e2_buf, int32_t* ptxt_buf,
                      double* r_tmp_d, double* e1_tmp_d, uint32_t* e2_tmp_u);

void Decaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* s_ptr, int32_t* z_ptr,
                      int32_t* c_in, std::vector<int32_t>& Hash_K,
                      const std::vector<int32_t>& h_c_in, /* NEW: reuse from Encaps */
                      int32_t* mp_buf, int32_t* dec_buf,
                      int32_t* r_buf, int32_t* e1_buf,
                      int32_t* e2_buf, int32_t* ptxt_buf, int32_t* c_chk,
                      double* r_tmp_d, double* e1_tmp_d, uint32_t* e2_tmp_u);

/* ---- on-the-fly A variants (unchanged signatures) ---- */
void Encaps_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* t_ptr, int32_t* c_out,
                std::vector<int32_t>& Hash_K,
                int32_t* m_buf, int32_t* r_buf, int32_t* e1_buf,
                int32_t* e2_buf, int32_t* ptxt_buf);

void Decaps_GPU(uint64_t key_seed, uint32_t n, uint32_t q,
                int32_t* t_ptr, int32_t* s_ptr, int32_t* z_ptr,
                int32_t* c_in, std::vector<int32_t>& Hash_K,
                int32_t* mp_buf, int32_t* dec_buf,
                int32_t* r_buf, int32_t* e1_buf,
                int32_t* e2_buf, int32_t* ptxt_buf, int32_t* c_chk);

/* ---- CPU reference ---- */
void Encaps(uint32_t n, uint32_t q,
            std::vector<int32_t>& t, std::vector<int32_t>& c,
            const std::vector<std::vector<int32_t>>& A,
            const std::vector<std::vector<int32_t>>& AT,
            std::vector<int32_t>& Hash_K);

void Decaps(uint32_t n, uint32_t q,
            const std::vector<int32_t>& t,
            const std::vector<int32_t>& s_k,
            const std::vector<int32_t>& c,
            std::vector<int32_t>& Hash_K,
            const std::vector<std::vector<int32_t>>& A,
            const std::vector<std::vector<int32_t>>& AT);
