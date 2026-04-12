#pragma once
#include <vector>
#include <cstdint>
#include <random>

#include <omp.h>

#include "rng_openssl.h"

// ============================================================================
// LEGACY MODULAR ARITHMETIC
// ============================================================================
int mod(int a, int b);

// ============================================================================
// OPTIMIZED FAST MODULAR ARITHMETIC
// ============================================================================
// Barrett reduction for any modulus
int32_t fast_mod(int64_t a, uint32_t q);

// Specialized Barrett reduction for q=3329 (most common case)
int32_t fast_mod_q3329(int64_t a);

// ============================================================================
// MATRIX AND VECTOR OPERATIONS
// ============================================================================
std::vector<std::vector<int32_t>> GenerateRandomMatrixInt32(size_t n, int32_t maxValue);

int32_t sample_eta_centered_binomial(uint8_t eta);

std::vector<int32_t> sample_vector_binomial(uint32_t n);

int32_t sample_discrete_gaussian(double sigma);

std::vector<int32_t> GenerateGaussianVector(size_t n);

int32_t getRandomInt(int min, int max);

int random_bit();

std::vector<int32_t> concat(const std::vector<int32_t> &a, const std::vector<int32_t> &b);

std::vector<int32_t> flatten_matrix(const std::vector<std::vector<int32_t>> &A);

void transposeA(const std::vector<std::vector<int32_t>> &A, std::vector<std::vector<int32_t>> &AT);

// ============================================================================
// OPTIMIZED OPERATIONS FOR PKE/KEM
// ============================================================================

// Fast dot product with SIMD acceleration
int64_t fast_dot_product(const std::vector<int32_t>& a, 
                         const std::vector<int32_t>& b,
                         size_t n);

// Optimized matrix-vector multiplication with parallelization
void matrix_vector_mul(const std::vector<std::vector<int32_t>>& A,
                       const std::vector<int32_t>& v,
                       std::vector<int32_t>& result,
                       size_t n,
                       uint32_t q);

