#pragma once
#include <vector>
#include <cstdint>
#include <random>

#include <omp.h> // openMP

// #include "rng.h"
#include "rng_openssl.h"
// #include "rng_keccak.h"

int mod(int a, int b);

std::vector<std::vector<int32_t>> GenerateRandomMatrixInt32(size_t n, int32_t maxValue);

int32_t sample_eta_centered_binomial(uint8_t eta, std::mt19937 &gen);

std::vector<int32_t> sample_vector_binomial(uint32_t n);

int32_t sample_discrete_gaussian(double sigma);

std::vector<int32_t> GenerateGaussianVector(size_t n);

int32_t getRandomInt(int min, int max);

int random_bit();

std::vector<int32_t> concat(const std::vector<int32_t> &a, const std::vector<int32_t> &b);

std::vector<int32_t> flatten_matrix(const std::vector<std::vector<int32_t>> &A);

void transposeA(const std::vector<std::vector<int32_t>> &A, std::vector<std::vector<int32_t>> &AT);