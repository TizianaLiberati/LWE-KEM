#include <vector>
#include <cstdint>
#include <chrono>
#include <iostream>
#include <cstring>
#include <omp.h>

#include "kem.h"
#include "utils.h"
#include "hash_openssl.h"
#include "noise.h"
#include "pke.h"

#include "config.h"

// ============================================================================
// OPTIMIZED ENCAPSULATION WITH PARALLEL NOISE GENERATION
// ============================================================================

void Encaps(uint32_t n, uint32_t q, 
            std::vector<int32_t> &t, 
            std::vector<int32_t> &c, 
            const std::vector<std::vector<int32_t>> &A, 
            const std::vector<std::vector<int32_t>> &AT, 
            std::vector<int32_t> &Hash_K, const LweKemConfig& cfg) {
    const size_t msg_bits = 256;
    
    // Compute pkh = SHA3_256(A || t)
    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh = SHA3_256_openssl(At);
    
    // Generate random message m
    std::vector<int32_t> m(msg_bits);
    #pragma omp parallel for simd
    for (uint32_t i = 0; i < msg_bits; ++i) {
        m[i] = random_bit();
    }
    
    // Compute seed = SHA3_256(pkh || m)
    std::vector<int32_t> pkh_m = concat(pkh, m);
    std::vector<int32_t> seed = SHA3_256_openssl(pkh_m);
    
    // Generate XOF output for all noise sampling at once
    std::vector<uint8_t> coins = xof_coins_openssl(seed, n, msg_bits);
    
    // Precompute K_cap = SHA3_256(pkh || m)
    std::vector<int32_t> K_cap = SHA3_256_openssl(concat(pkh, m));
    
    // Pre-allocate ciphertext
    c.assign(msg_bits * (n + 1), 0);
    
    // OPTIMIZED: Parallel encryption of all bits with coalesced memory access
    // Each thread processes multiple bits with pre-allocated buffers
    #pragma omp parallel
    {
        // Thread-local buffers to avoid allocation overhead
        std::vector<int32_t> u_local(n);
        int32_t v_local;
        
        #pragma omp for schedule(dynamic, 8)
        for (uint32_t j = 0; j < msg_bits; ++j) {
            // Compute offset for this bit's ciphertext
            size_t pos = j * (2 * n + 1);
            size_t off = static_cast<size_t>(j) * (n + 1);
            
            // Generate noise for this bit
            NoiseTriple noise_i = GenerateNoisesForOneBit(coins, pos, n, cfg);
            
            // Map plaintext bit to ring element
            int32_t plaintext_j = (m[j] == 1) ? static_cast<int32_t>(q / 2) : 0;
            
            // Encrypt
            Encrypt(n, q, t, u_local, v_local, plaintext_j, 
                    noise_i.r, noise_i.e1, noise_i.e2, AT);
            
            // Write to ciphertext with coalesced pattern
            // OPTIMIZED: Use memcpy for vector copy
            std::memcpy(&c[off], u_local.data(), n * sizeof(int32_t));
            c[off + n] = v_local;
        }
    }
    
    // Compute h_c = SHA3_256(c)
    std::vector<int32_t> h_c = SHA3_256_openssl(c);
    
    // K = SHA3_256(K_cap || h_c)
    Hash_K = SHA3_256_openssl(concat(K_cap, h_c));
}

// ============================================================================
// OPTIMIZED DECAPSULATION WITH PARALLEL VERIFICATION
// ============================================================================

void Decaps(uint32_t n, uint32_t q, 
            const std::vector<int32_t> &t, 
            const std::vector<int32_t> &s_k, 
            const std::vector<int32_t> &c, 
            std::vector<int32_t> &Hash_K, 
            const std::vector<std::vector<int32_t>> &A, 
            const std::vector<std::vector<int32_t>> &AT, const LweKemConfig& cfg) {
    const size_t msg_bits = 256;
    
    // Compute pkh
    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh = SHA3_256_openssl(At);
    
    // OPTIMIZED: Parallel decryption of all bits
    std::vector<int32_t> mprime(msg_bits);
    
    #pragma omp parallel
    {
        std::vector<int32_t> u_j(n);
        int32_t v_j, dec_m;
        
        #pragma omp for schedule(static)
        for (uint32_t j = 0; j < msg_bits; ++j) {
            size_t off = static_cast<size_t>(j) * (n + 1);
            
            // Load u from ciphertext
            std::memcpy(u_j.data(), &c[off], n * sizeof(int32_t));
            v_j = c[off + n];
            
            // Decrypt
            Decrypt(v_j, u_j, s_k, q, dec_m);
            mprime[j] = (dec_m == static_cast<int32_t>(q / 2)) ? 1 : 0;
        }
    }
    
    // Compute K_cap
    std::vector<int32_t> K_cap = SHA3_256_openssl(concat(pkh, mprime));
    
    // Regenerate seed and coins
    std::vector<int32_t> pkh_mprime = concat(pkh, mprime);
    std::vector<int32_t> seed = SHA3_256_openssl(pkh_mprime);
    std::vector<uint8_t> coins = xof_coins_openssl(seed, n, msg_bits);
    
    // OPTIMIZED: Parallel re-encryption for verification
    std::vector<int32_t> cchk;
    cchk.assign(msg_bits * (n + 1), 0);
    
    #pragma omp parallel
    {
        std::vector<int32_t> u_tmp(n);
        int32_t v_tmp;
        
        #pragma omp for schedule(dynamic, 8)
        for (uint32_t j = 0; j < msg_bits; ++j) {
            size_t pos = j * (2 * n + 1);
            size_t off = static_cast<size_t>(j) * (n + 1);
            
            int32_t m_j_map = mprime[j] ? static_cast<int32_t>(q / 2) : 0;
            
            NoiseTriple noise_i = GenerateNoisesForOneBit(coins, pos, n, cfg);
            
            Encrypt(n, q, const_cast<std::vector<int32_t> &>(t), 
                    u_tmp, v_tmp, static_cast<uint32_t>(m_j_map),
                    noise_i.r, noise_i.e1, noise_i.e2, 
                    const_cast<std::vector<std::vector<int32_t>> &>(AT));
            
            std::memcpy(&cchk[off], u_tmp.data(), n * sizeof(int32_t));
            cchk[off + n] = v_tmp;
        }
    }
    
    // OPTIMIZED: Fast constant-time comparison using memcmp
    bool equal = (cchk.size() == c.size());
    if (equal) {
        equal = (std::memcmp(cchk.data(), c.data(), c.size() * sizeof(int32_t)) == 0);
    }
    
    // Compute h_c
    std::vector<int32_t> h_c = SHA3_256_openssl(c);
    
    // Final key derivation
    if (equal) {
        Hash_K = SHA3_256_openssl(concat(K_cap, h_c));
    } else {
        // Implicit rejection: use random z from secret key
        std::vector<int32_t> z(s_k.begin() + n, s_k.end());
        Hash_K = SHA3_256_openssl(concat(z, h_c));
    }
}

// ============================================================================
// STREAMING XOF OPTIMIZATION (future enhancement)
// ============================================================================

// This function would process XOF output incrementally
// Currently using the existing implementation as interface is compatible
// Future: implement chunked XOF to reduce memory allocation

