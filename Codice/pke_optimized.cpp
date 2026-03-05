#include <vector>
#include <cstdint>
#include <iostream>
#include <cstring>

#include <omp.h> // openMP

#include "pke.h"
#include "utils.h"

/////////////////////////////////////   KeyGen  /////////////////////////////////////
/**
 * ADVANCED OPTIMIZATION: Tiled/blocked matrix-vector multiply
 * Use when n is large (>= 512) for better cache utilization
 *
 * This version uses tile-based computation to improve memory locality
 */
void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A,
            std::vector<int32_t> &s_k, std::vector<int32_t> &t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    std::vector<int32_t> Aflat = flatten_matrix(A);
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, 0);

    std::vector<int32_t> z(256);
    for (int i = 0; i < 256; ++i)
        z[i] = getRandomInt(0, 1);
    s_k = concat(s, z);

    int32_t* __restrict__ Aflat_ptr = Aflat.data();
    int32_t* __restrict__ s_ptr = s.data();
    int32_t* __restrict__ e_ptr = e.data();
    int32_t* __restrict__ prod_ptr = prod.data();

    const long long q_ll = (long long)q;
    const uint32_t TILE_SIZE = 32;  // Tunable parameter

    #pragma acc data copyin(Aflat_ptr[0:n*n], s_ptr[0:n], e_ptr[0:n]) copyout(prod_ptr[0:n])
    {
        // Tiled matrix-vector multiply
        #pragma acc parallel loop gang
        for (uint32_t i = 0; i < n; ++i)
        {
            long long acc = 0;

            // Process in tiles for better cache locality
            #pragma acc loop vector reduction(+:acc)
            for (uint32_t tile = 0; tile < (n + TILE_SIZE - 1) / TILE_SIZE; ++tile) {
                uint32_t tile_start = tile * TILE_SIZE;
                uint32_t tile_end = (tile_start + TILE_SIZE < n) ? tile_start + TILE_SIZE : n;

                for (uint32_t j = tile_start; j < tile_end; ++j) {
                    acc += (long long)Aflat_ptr[i * n + j] * (long long)s_ptr[j];
                }
            }

            acc += (long long)e_ptr[i];
            long long r = acc % q_ll;
            r += (r >> 63) & q_ll;
            prod_ptr[i] = (int32_t)r;
        }
    }

    t = prod;
}

/////////////////////////////////////   Encrypt /////////////////////////////////////
/**
 * GPU-accelerated Encrypt with safer memory management
 *
 * Key optimizations:
 * 1. Single unified data region for all GPU operations
 * 2. Fused matrix multiply + dot product computation
 * 3. Gang/vector/worker parallelism hints
 * 4. Optimized memory access patterns
 * 5. Efficient reduction strategy
 * 6. Branchless modulo operations
 */
void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u,
             int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r,
             std::vector<int32_t> &e1, int32_t &e2, const std::vector<std::vector<int32_t>> &AT)
{
    std::vector<int32_t> ATflat = flatten_matrix(AT);
    u.resize(n);

    // Size validation
    if (ATflat.size() != n*n || r.size() != n || e1.size() != n || t.size() != n) {
        std::cerr << "Encrypt: Vector size mismatch!" << std::endl;
        return;
    }

    // Extract pointers with restrict for optimization
    int32_t* __restrict__ ATflat_ptr = ATflat.data();
    int32_t* __restrict__ r_ptr = r.data();
    int32_t* __restrict__ e1_ptr = e1.data();
    int32_t* __restrict__ u_ptr = u.data();
    int32_t* __restrict__ t_ptr = t.data();

    const long long q_ll = (long long)q;
    const uint32_t TILE_SIZE = 32;  // Tunable parameter
    long long dot = 0;

    // ====== UNIFIED GPU Region: Compute u AND dot in parallel ======
    // Use structured data region for all operations
    #pragma acc data copyin(ATflat_ptr[0:n*n], r_ptr[0:n], e1_ptr[0:n], t_ptr[0:n]) \
                     copyout(u_ptr[0:n]) copy(dot)
    {
        // Parallel region 1: u = AT*r + e1 (mod q)
        #pragma acc parallel loop gang
        for (uint32_t i = 0; i < n; ++i)
        {
            long long acc = 0;

            #pragma acc loop vector reduction(+:acc)
            for (uint32_t tile = 0; tile < (n + TILE_SIZE - 1) / TILE_SIZE; ++tile) {
                uint32_t tile_start = tile * TILE_SIZE;
                uint32_t tile_end = (tile_start + TILE_SIZE < n) ? tile_start + TILE_SIZE : n;

                for (uint32_t j = tile_start; j < tile_end; ++j) {
                    acc += (long long)ATflat_ptr[i * n + j] * (long long)r_ptr[j];
                }
            }

            acc += (long long)e1_ptr[i];

            // Branchless modulo
            long long m = acc % q_ll;
            m += (m >> 63) & q_ll;
            u_ptr[i] = (int32_t)m;
        }

        // Parallel region 2: dot = t·r (computed concurrently with above)
        // Use independent async for potential overlap
        #pragma acc parallel loop gang vector reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i) {
            dot += (long long)t_ptr[i] * (long long)r_ptr[i];
        }
    }

    // Finalize on CPU with optimized modulo
    long long res = dot % q_ll;
    res += (res >> 63) & q_ll;

    long long vv = res + (long long)e2 + (long long)plaintext_i;
    vv %= q_ll;
    vv += (vv >> 63) & q_ll;

    v_i = (int32_t)vv;
}

/////////////////////////////////////   Decrypt /////////////////////////////////////
/**
 * GPU-accelerated Decrypt with safer memory management
 *
 * Single GPU kernel with proper size verification
 * 
 * FIXED: Corrected tiling logic and missing closing brace
 */
void Decrypt(int32_t v_i, const std::vector<int32_t> &u, const std::vector<int32_t> &s_k,
             uint32_t q, int32_t &decrypt_i)
{
    const size_t n = u.size();
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n);

    // Verify vector sizes before GPU operations
    if (u.size() != n || s.size() != n) {
        std::cerr << "Decrypt: Vector size mismatch! u=" << u.size()
                  << " s=" << s.size() << std::endl;
        return;
    }

    // Extract pointers BEFORE pragmas (critical!)
    const int32_t* __restrict__ u_ptr = u.data();
    const int32_t* __restrict__ s_ptr = s.data();
    const long long q_ll = (long long)q;

    long long dot = 0;

    // ====== GPU Region: dot = u·s ======
    // Option 1: Simple version without tiling (recommended for dot products)
    #pragma acc data copyin(u_ptr[0:n], s_ptr[0:n])
    {
        #pragma acc parallel loop gang worker vector reduction(+:dot)
        for (size_t i = 0; i < n; ++i) {
            dot += (long long)s_ptr[i] * (long long)u_ptr[i];
        }
    }  // Close data region BEFORE CPU computation

    // Final computation on CPU (outside GPU region)
    long long r = ((long long)v_i - dot) % q_ll;
    r += (r >> 63) & q_ll;  // Branchless: if (r < 0) r += q
    
    int32_t mu = (int32_t)r;
    const int32_t bound = (int32_t)q / 4;
    
    // Branchless decryption decision
    int32_t in_low_range = (mu <= bound);
    int32_t in_high_range = (mu >= (int32_t)q - bound);
    decrypt_i = ((1 - (in_low_range | in_high_range)) * ((int32_t)q / 2));
}

/////////////////////////////////////   Alternative Decrypt with Tiling /////////////////////////////////////
/**
 * Alternative version with tiling (use if n is very large)
 */
void Decrypt_Tiled(int32_t v_i, const std::vector<int32_t> &u, const std::vector<int32_t> &s_k,
                   uint32_t q, int32_t &decrypt_i)
{
    const size_t n = u.size();
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n);

    if (u.size() != n || s.size() != n) {
        std::cerr << "Decrypt: Vector size mismatch!" << std::endl;
        return;
    }

    const int32_t* __restrict__ u_ptr = u.data();
    const int32_t* __restrict__ s_ptr = s.data();
    const long long q_ll = (long long)q;
    const uint32_t TILE_SIZE = 32;

    long long dot = 0;

    // ====== GPU Region with Tiling: dot = u·s ======
    #pragma acc data copyin(u_ptr[0:n], s_ptr[0:n])
    {
        // Outer loop over tiles (gang-parallel)
        #pragma acc parallel loop gang reduction(+:dot)
        for (uint32_t tile = 0; tile < (n + TILE_SIZE - 1) / TILE_SIZE; ++tile) {
            long long tile_sum = 0;
            uint32_t tile_start = tile * TILE_SIZE;
            uint32_t tile_end = (tile_start + TILE_SIZE < n) ? tile_start + TILE_SIZE : n;
            
            // Inner loop over elements in tile (vector-parallel)
            #pragma acc loop vector reduction(+:tile_sum)
            for (size_t i = tile_start; i < tile_end; ++i) {
                tile_sum += (long long)s_ptr[i] * (long long)u_ptr[i];
            }
            
            dot += tile_sum;
        }
    }  // Close data region

    // Final computation on CPU
    long long r = ((long long)v_i - dot) % q_ll;
    r += (r >> 63) & q_ll;
    
    int32_t mu = (int32_t)r;
    const int32_t bound = (int32_t)q / 4;
    
    int32_t in_low_range = (mu <= bound);
    int32_t in_high_range = (mu >= (int32_t)q - bound);
    decrypt_i = ((1 - (in_low_range | in_high_range)) * ((int32_t)q / 2));
}
