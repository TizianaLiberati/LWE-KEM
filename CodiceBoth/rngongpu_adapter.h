#pragma once

#include <cstdint>
#include <cstddef>   // for size_t
#include <vector>

/* ============================================================================
 *  GPU / RNGonGPU adapter interface (C++ linkage)
 *
 *  IMPORTANT:
 *  - These functions are implemented across .cpp / .cu files compiled as C++
 *  - DO NOT use extern "C" here (would break name mangling)
 * ============================================================================ */

/* ===== Matrix generation ===== */
void GenerateA_GPU_rngongpu(uint64_t rho_seed,
                           uint32_t n,
                           uint32_t q,
                           int32_t* Ap);

/* ===== Key generation ===== */
void KeyGen_GPU_rngongpu_Aflat(uint64_t key_seed,
                              uint32_t n,
                              uint32_t q,
                              int32_t* Ap,
                              int32_t* sp,
                              int32_t* tp,
                              uint32_t* sraw_p,
                              double* etmp_p,
                              int32_t* ebuf_p);

/* ===== Noise generation (GPU) ===== */
void GPU_GenerateNoise_rngongpu(uint64_t noise_seed,
                               uint32_t n,
                               uint32_t msg_bits,
                               int32_t* rp,
                               int32_t* e1p,
                               int32_t* e2p,
                               double* r_tmp,
                               double* e1_tmp,
                               uint32_t* e2_tmp);

/* ===== Batch encryption (GPU) ===== */
void BatchEncrypt_GPU_Aflat(uint32_t n,
                           uint32_t q,
                           const int32_t* A_flat,
                           const int32_t* t_ptr,
                           const int32_t* r,
                           const int32_t* e1,
                           const int32_t* e2,
                           const int32_t* m,
                           int32_t* c_out,
                           uint32_t msg_bits);

/* ===== Ciphertext comparison ===== */
bool GPU_CiphertextEqual(const int32_t* c1,
                         const int32_t* c2,
                         size_t len);

/* ===== Encapsulation ===== */
void Encaps_GPU_Aflat(uint64_t key_seed,
                      uint32_t n,
                      uint32_t q,
                      int32_t* Ap,
                      int32_t* tp,
                      int32_t* cp,
                      std::vector<int32_t>& K_enc,
                      std::vector<int32_t>& h_c_enc,
                      int32_t* mp,
                      int32_t* rp,
                      int32_t* e1p,
                      int32_t* e2p,
                      int32_t* ptxtp,
                      double* r_tmp_d,
                      double* e1_tmp_d,
                      uint32_t* e2_tmp_u);

/* ===== Decapsulation ===== */
void Decaps_GPU_Aflat(uint64_t key_seed,
                      uint32_t n,
                      uint32_t q,
                      int32_t* Ap,
                      int32_t* tp,
                      int32_t* sp,
                      int32_t* zp,
                      int32_t* cp,
                      std::vector<int32_t>& K_dec,
                      const std::vector<int32_t>& h_c_enc,
                      int32_t* mpp,
                      int32_t* decp,
                      int32_t* rp,
                      int32_t* e1p,
                      int32_t* e2p,
                      int32_t* ptxtp,
                      int32_t* cchkp,
                      double* r_tmp_d,
                      double* e1_tmp_d,
                      uint32_t* e2_tmp_u);
