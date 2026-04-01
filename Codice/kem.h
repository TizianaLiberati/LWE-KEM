#pragma once
#include <vector>
#include <cstdint>
#include "utils.h"
#include "hash_openssl.h"
#include "noise.h"
#include "pke.h"
#include "rng_openssl.h"

/* Original CPU versions */
void Encaps(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &c,
            const std::vector<std::vector<int32_t>> &A,
            const std::vector<std::vector<int32_t>> &AT,
            std::vector<int32_t> &Hash_K);

void Decaps(uint32_t n, uint32_t q, const std::vector<int32_t> &t,
            const std::vector<int32_t> &s_k, const std::vector<int32_t> &c,
            std::vector<int32_t> &Hash_K,
            const std::vector<std::vector<int32_t>> &A,
            const std::vector<std::vector<int32_t>> &AT);

/* RNGonGPU verions*/
void Encaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* c_out,
                      std::vector<int32_t>& Hash_K,
                      int32_t* m_buf, int32_t* r_buf, int32_t* e1_buf,
                      int32_t* e2_buf, int32_t* ptxt_buf);

void Decaps_GPU_Aflat(uint64_t key_seed, uint32_t n, uint32_t q,
                      int32_t* A_flat, int32_t* t_ptr, int32_t* s_ptr, int32_t* z_ptr,
                      int32_t* c_in, std::vector<int32_t>& Hash_K,
                      int32_t* mp_buf, int32_t* dec_buf,
                      int32_t* r_buf, int32_t* e1_buf,
                      int32_t* e2_buf, int32_t* ptxt_buf, int32_t* c_chk);            

/* GPU-accelerated versions */
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