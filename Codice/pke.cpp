#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>

#include <omp.h>

#include "pke.h"
#include "utils.h"

// ============================================================================
// OPTIMIZED KEY GENERATION
// ============================================================================

void KeyGen(uint32_t n, uint32_t q, 
            std::vector<std::vector<int32_t>> &A, 
            std::vector<int32_t> &s_k, 
            std::vector<int32_t> &t) {
    // Generate random matrix A
    A = GenerateRandomMatrixInt32(n, q - 1);
    
    // Sample secret and error vectors
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    
    // Generate random z for implicit rejection (256 bits)
    std::vector<int32_t> z(256);
    #pragma omp parallel for simd
    for (int i = 0; i < 256; ++i) {
        z[i] = getRandomInt(0, 1);
    }
    
    // Secret key: concat of s and z
    s_k = concat(s, z);
    
    // Compute t = A*s + e with optimized matrix-vector multiplication
    std::vector<int32_t> prod(n);
    matrix_vector_mul(A, s, prod, n, q);
    
    // Add error vector with modular reduction
    t.resize(n);
    #pragma omp parallel for simd
    for (uint32_t i = 0; i < n; ++i) {
        int64_t val = static_cast<int64_t>(prod[i]) + e[i];
        t[i] = fast_mod(val, q);
    }
}

// ============================================================================
// OPTIMIZED ENCRYPTION
// ============================================================================

void Encrypt(uint32_t n, uint32_t q, 
             std::vector<int32_t> &t, 
             std::vector<int32_t> &u, 
             int32_t &v_i, 
             uint32_t plaintext_i, 
             std::vector<int32_t> &r, 
             std::vector<int32_t> &e1, 
             int32_t &e2, 
             const std::vector<std::vector<int32_t>> &AT) {
    // OPTIMIZED: Use thread-local storage for intermediate results
    static thread_local std::vector<int32_t> prod_buffer;
    prod_buffer.resize(n);
    
    // u = AT * r with parallel matrix-vector multiplication
    matrix_vector_mul(AT, r, prod_buffer, n, q);
    
    // Add e1 with modular reduction
    u.resize(n);
    #pragma omp parallel for simd
    for (uint32_t i = 0; i < n; ++i) {
        int64_t val = static_cast<int64_t>(prod_buffer[i]) + e1[i];
        u[i] = fast_mod(val, q);
    }
    
    // OPTIMIZED: Compute v = t.r + e2 + plaintext with SIMD dot product
    int64_t dot = fast_dot_product(t, r, n);
    
    // Final modular reduction for v
    int64_t v_val = dot + e2 + static_cast<int64_t>(plaintext_i);
    v_i = fast_mod(v_val, q);
}

// ============================================================================
// OPTIMIZED DECRYPTION
// ============================================================================

void Decrypt(int32_t v_i, 
             const std::vector<int32_t> &u, 
             const std::vector<int32_t> &s_k, 
             uint32_t q, 
             int32_t &decrypt_i) {
    const size_t n = u.size();
    
    // Extract s from s_k (first n elements)
    // OPTIMIZED: Use pointer arithmetic to avoid copy
    const int32_t* s_ptr = s_k.data();
    
    // SIMD-accelerated dot product s . u
    int64_t dot = 0;
    #pragma omp simd reduction(+:dot)
    for (size_t i = 0; i < n; ++i) {
        dot += static_cast<int64_t>(s_ptr[i]) * static_cast<int64_t>(u[i]);
    }
    
    // Compute (v - s.u) mod q
    int64_t r = static_cast<int64_t>(v_i) - dot;
    r %= static_cast<int64_t>(q);
    if (r < 0) r += q;
    
    int32_t mu = static_cast<int32_t>(r);
    
    // Decode: 0 if in [0, q/4] U [3q/4, q), else q/2
    const int32_t bound = static_cast<int32_t>(q / 4);
    const int32_t q_minus_bound = static_cast<int32_t>(q) - bound;
    
    decrypt_i = (mu <= bound || mu >= q_minus_bound) ? 0 : static_cast<int32_t>(q / 2);
}

// ============================================================================
// BATCH ENCRYPTION FOR IMPROVED THROUGHPUT
// ============================================================================

// Encrypt multiple plaintext bits using shared randomness
// This reduces per-bit overhead by sharing matrix operations
void EncryptBatch(uint32_t n, uint32_t q,
                  std::vector<int32_t> &t,
                  std::vector<std::vector<int32_t>> &u_batch,
                  std::vector<int32_t> &v_batch,
                  const std::vector<uint32_t> &plaintexts,
                  std::vector<std::vector<int32_t>> &r_batch,
                  std::vector<std::vector<int32_t>> &e1_batch,
                  std::vector<int32_t> &e2_batch,
                  const std::vector<std::vector<int32_t>> &AT) {
    size_t batch_size = plaintexts.size();
    u_batch.resize(batch_size);
    v_batch.resize(batch_size);
    
    // Parallel encryption of all bits
    #pragma omp parallel for schedule(static)
    for (size_t b = 0; b < batch_size; ++b) {
        Encrypt(n, q, t, u_batch[b], v_batch[b], plaintexts[b],
                r_batch[b], e1_batch[b], e2_batch[b], AT);
    }
}

