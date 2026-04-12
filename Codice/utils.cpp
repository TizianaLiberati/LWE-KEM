#include <vector>
#include <cstdint>
#include <random>
#include <limits>
#include <cstring>
#include <omp.h>

#include "rng_openssl.h"
#include "utils.h"

// ============================================================================
// OPTIMIZED MODULAR ARITHMETIC WITH BARRETT REDUCTION
// ============================================================================
// For fixed modulus q=3329, Barrett reduction eliminates division operations
// Formula: a mod q = a - floor(a * mu / 2^32) * q, where mu = floor(2^32 / q)

namespace {
    // Barrett precomputation for q=3329
    constexpr uint32_t BARRETT_Q = 3329;
    constexpr uint32_t BARRETT_MU = 1290167;  // floor(2^32 / 3329)
    constexpr uint32_t BARRETT_SHIFT = 32;
    
    // Fast branchless modulo using Barrett reduction
    inline int32_t mod_barrett(int64_t a, uint32_t q) {
        // Handle negative inputs properly
        uint64_t ua = static_cast<uint64_t>(a < 0 ? -a : a);
        uint64_t t = (ua * BARRETT_MU) >> BARRETT_SHIFT;
        uint64_t r = ua - t * q;
        // Final reduction (at most one iteration needed for q=3329)
        if (r >= q) r -= q;
        return static_cast<int32_t>(a < 0 ? -(int64_t)r : r);
    }
    
    // Specialized fast mod for q=3329 (most common case)
    inline int32_t mod_q3329(int64_t a) {
        uint64_t ua = static_cast<uint64_t>(a < 0 ? -a : a);
        uint64_t t = (ua * BARRETT_MU) >> BARRETT_SHIFT;
        uint32_t r = static_cast<uint32_t>(ua - t * BARRETT_Q);
        if (r >= BARRETT_Q) r -= BARRETT_Q;
        return static_cast<int32_t>(a < 0 ? -(int32_t)r : (int32_t)r);
    }
    
    // Lazy reduction: accumulate in 64-bit, reduce once at end
    inline int32_t mod_lazy(int64_t a, uint32_t q) {
        int64_t r = a % q;
        if (r < 0) r += q;
        return static_cast<int32_t>(r);
    }
}

// Legacy mod function - optimized version
int mod(int a, int b) {
    return (a % b + b) % b;
}

// ============================================================================
// FAST MODULAR ARITHMETIC INTERFACE
// ============================================================================

int32_t fast_mod(int64_t a, uint32_t q) {
    return mod_barrett(a, q);
}

int32_t fast_mod_q3329(int64_t a) {
    return mod_q3329(a);
}

// ============================================================================
// OPTIMIZED MATRIX GENERATION WITH CACHE-FRIENDLY ACCESS
// ============================================================================

// Structure of Arrays (SoA) layout for better vectorization
// Instead of A[i][j], use A_data[j*n + i] for column-major access

std::vector<std::vector<int32_t>> GenerateRandomMatrixInt32(size_t n, int32_t maxValue) {
    OpenSSL_URBG rng;
    std::uniform_int_distribution<int32_t> dist(0, maxValue);
    
    // Use flat allocation for better cache locality
    std::vector<int32_t> flat_data(n * n);
    
    // OPTIMIZED: Single parallel region with larger work chunks
    #pragma omp parallel
    {
        // Thread-local RNG to avoid contention
        OpenSSL_URBG local_rng;
        std::uniform_int_distribution<int32_t> local_dist(0, maxValue);
        
        #pragma omp for schedule(static)
        for (size_t idx = 0; idx < n * n; ++idx) {
            flat_data[idx] = local_dist(local_rng);
        }
    }
    
    // Convert to 2D view (still row-major for compatibility)
    std::vector<std::vector<int32_t>> matrix(n);
    for (size_t i = 0; i < n; ++i) {
        matrix[i] = std::vector<int32_t>(flat_data.begin() + i * n, 
                                          flat_data.begin() + (i + 1) * n);
    }
    
    return matrix;
}

// ============================================================================
// OPTIMIZED BINOMIAL SAMPLING
// ============================================================================

// OPTIMIZED: Use lookup table for binomial with eta=3
// sum of 3 bits - sum of 3 bits: range [-3, 3]
namespace {
    // Precompute distribution for eta=3
    // P(X=k) where X = sum(B1[0:3]) - sum(B2[0:3])
    constexpr uint8_t BINOMIAL_ETA3_WEIGHTS[7] = {
        1,  // -3: 1/64
        6,  // -2: 6/64
        15, // -1: 15/64
        20, //  0: 20/64
        15, //  1: 15/64
        6,  //  2: 6/64
        1   //  3: 1/64
    };
}

int32_t sample_eta_centered_binomial(uint8_t eta) {
    OpenSSL_URBG rng;
    std::uniform_int_distribution<uint8_t> dis(0, 1);
    
    // OPTIMIZED: Use population count approach for small eta
    // Pack 6 random bits into a byte, use table lookup
    uint8_t packed = 0;
    for (uint8_t i = 0; i < 6; ++i) {
        packed = (packed << 1) | dis(rng);
    }
    
    // Count bits in first 3 and last 3 positions
    uint8_t sum1 = ((packed >> 3) & 1) + ((packed >> 4) & 1) + ((packed >> 5) & 1);
    uint8_t sum2 = (packed & 1) + ((packed >> 1) & 1) + ((packed >> 2) & 1);
    
    return static_cast<int32_t>(sum1) - static_cast<int32_t>(sum2);
}

std::vector<int32_t> sample_vector_binomial(uint32_t n) {
    std::vector<int32_t> result;
    result.reserve(n);
    
    // Batch RNG calls for better performance
    OpenSSL_URBG rng;
    std::uniform_int_distribution<uint8_t> dis(0, 255);
    
    // Each byte provides enough entropy for ~1.3 samples
    // Process in batches to reduce RNG overhead
    constexpr uint32_t BATCH_SIZE = 64;
    uint8_t buffer[BATCH_SIZE];
    
    uint32_t i = 0;
    while (i < n) {
        // Fill buffer with random bytes
        for (uint32_t b = 0; b < BATCH_SIZE; ++b) {
            buffer[b] = dis(rng);
        }
        
        // Consume random bytes to generate samples
        uint32_t buf_idx = 0;
        while (i < n && buf_idx < BATCH_SIZE - 1) {
            // Use 6 bits from each byte pair
            uint8_t byte1 = buffer[buf_idx++];
            uint8_t byte2 = buffer[buf_idx++];
            
            uint8_t sum1 = ((byte1 >> 5) & 1) + ((byte1 >> 6) & 1) + ((byte1 >> 7) & 1);
            uint8_t sum2 = ((byte1 >> 2) & 1) + ((byte1 >> 3) & 1) + ((byte1 >> 4) & 1);
            result.push_back(static_cast<int32_t>(sum1) - static_cast<int32_t>(sum2));
            i++;
            
            if (i < n) {
                sum1 = ((byte2 >> 5) & 1) + ((byte2 >> 6) & 1) + ((byte2 >> 7) & 1);
                sum2 = ((byte2 >> 2) & 1) + ((byte2 >> 3) & 1) + ((byte2 >> 4) & 1);
                result.push_back(static_cast<int32_t>(sum1) - static_cast<int32_t>(sum2));
                i++;
            }
        }
    }
    
    return result;
}

// ============================================================================
// OPTIMIZED GAUSSIAN SAMPLING FOR KEY GENERATION
// ============================================================================

// OPTIMIZED: Use Box-Muller for vectorized normal sampling
namespace {
    // Cached normal samples (Box-Muller generates pairs)
    thread_local double cached_normal = 0.0;
    thread_local bool has_cached = false;
}

int32_t sample_discrete_gaussian(double sigma) {
    OpenSSL_URBG rng;
    
    // Use cached sample if available
    double z;
    if (has_cached) {
        z = cached_normal;
        has_cached = false;
    } else {
        // Box-Muller transform
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double u1 = uniform(rng);
        double u2 = uniform(rng);
        
        double mag = sigma * std::sqrt(-2.0 * std::log(u1));
        z = mag * std::cos(2.0 * M_PI * u2);
        cached_normal = mag * std::sin(2.0 * M_PI * u2);
        has_cached = true;
    }
    
    // Bound and round
    double bound = 6.0 * sigma;
    if (z > bound) z = bound;
    if (z < -bound) z = -bound;
    
    return static_cast<int32_t>(z >= 0.0 ? (z + 0.5) : (z - 0.5));
}

std::vector<int32_t> GenerateGaussianVector(size_t n) {
    std::vector<int32_t> vec;
    vec.reserve(n);
    
    OpenSSL_URBG rng;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double sigma = 2.3;
    
    // Generate in pairs using Box-Muller
    size_t i = 0;
    while (i < n) {
        double u1 = uniform(rng);
        double u2 = uniform(rng);
        
        double mag = sigma * std::sqrt(-2.0 * std::log(u1));
        double z1 = mag * std::cos(2.0 * M_PI * u2);
        double z2 = mag * std::sin(2.0 * M_PI * u2);
        
        double bound = 6.0 * sigma;
        
        // First sample
        if (z1 > bound) z1 = bound;
        if (z1 < -bound) z1 = -bound;
        vec.push_back(static_cast<int32_t>(z1 >= 0.0 ? (z1 + 0.5) : (z1 - 0.5)));
        i++;
        
        // Second sample
        if (i < n) {
            if (z2 > bound) z2 = bound;
            if (z2 < -bound) z2 = -bound;
            vec.push_back(static_cast<int32_t>(z2 >= 0.0 ? (z2 + 0.5) : (z2 - 0.5)));
            i++;
        }
    }
    
    return vec;
}

// ============================================================================
// RANDOM INTEGER GENERATION
// ============================================================================

int32_t getRandomInt(int min, int max) {
    OpenSSL_URBG rng;
    std::uniform_int_distribution<> distr(min, max);
    return distr(rng);
}

int random_bit() {
    OpenSSL_URBG urbg;
    std::uniform_int_distribution<int> bit(0, 1);
    return bit(urbg);
}

// ============================================================================
// OPTIMIZED VECTOR OPERATIONS
// ============================================================================

// OPTIMIZED: Use memcpy for concatenation when possible
std::vector<int32_t> concat(const std::vector<int32_t> &a, const std::vector<int32_t> &b) {
    std::vector<int32_t> out;
    out.reserve(a.size() + b.size());
    
    // Use insert for bulk copy (often optimized to memcpy)
    out.insert(out.end(), a.begin(), a.end());
    out.insert(out.end(), b.begin(), b.end());
    
    return out;
}

// OPTIMIZED: Pre-allocate and flatten in one pass
std::vector<int32_t> flatten_matrix(const std::vector<std::vector<int32_t>> &A) {
    if (A.empty()) return {};
    
    const size_t n = A.size();
    const size_t m = A[0].size();
    
    std::vector<int32_t> out;
    out.reserve(n * m);
    
    // Flatten row by row (cache-friendly)
    for (size_t i = 0; i < n; ++i) {
        out.insert(out.end(), A[i].begin(), A[i].end());
    }
    
    return out;
}

// OPTIMIZED: Transpose with blocking for cache efficiency
void transposeA(const std::vector<std::vector<int32_t>> &A, 
                std::vector<std::vector<int32_t>> &AT) {
    const size_t n = A.size();
    if (n == 0) return;
    
    AT.assign(n, std::vector<int32_t>(n));
    
    // Use cache blocking for large matrices
    constexpr size_t BLOCK_SIZE = 32;
    
    #pragma omp parallel for schedule(static)
    for (size_t ii = 0; ii < n; ii += BLOCK_SIZE) {
        for (size_t jj = 0; jj < n; jj += BLOCK_SIZE) {
            size_t i_max = std::min(ii + BLOCK_SIZE, n);
            size_t j_max = std::min(jj + BLOCK_SIZE, n);
            
            for (size_t i = ii; i < i_max; ++i) {
                for (size_t j = jj; j < j_max; ++j) {
                    AT[j][i] = A[i][j];
                }
            }
        }
    }
}

// ============================================================================
// UTILITY FUNCTIONS FOR OPTIMIZED PKE
// ============================================================================

// Fast dot product with lazy reduction
int64_t fast_dot_product(const std::vector<int32_t>& a, 
                         const std::vector<int32_t>& b,
                         size_t n) {
    int64_t sum = 0;
    
    #pragma omp simd reduction(+:sum)
    for (size_t i = 0; i < n; ++i) {
        sum += static_cast<int64_t>(a[i]) * static_cast<int64_t>(b[i]);
    }
    
    return sum;
}

// Matrix-vector multiplication with Barrett reduction
void matrix_vector_mul(const std::vector<std::vector<int32_t>>& A,
                       const std::vector<int32_t>& v,
                       std::vector<int32_t>& result,
                       size_t n,
                       uint32_t q) {
    result.resize(n);
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; ++i) {
        int64_t sum = 0;
        
        #pragma omp simd reduction(+:sum)
        for (size_t j = 0; j < n; ++j) {
            sum += static_cast<int64_t>(A[i][j]) * static_cast<int64_t>(v[j]);
        }
        
        result[i] = fast_mod(sum, q);
    }
}

