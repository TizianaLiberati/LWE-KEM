/*  utils.cpp  –  CPU-side utility functions
 *
 *  All functions run on CPU.  With -gpu=managed, the vectors they
 *  produce live in unified memory and are automatically accessible
 *  from GPU kernels.
 */

#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>

#include "rng_openssl.h"
#include "utils.h"

int mod(int a, int b)
{
    return (a % b + b) % b;
}

/////////////////////////////////////   Matrix / transpose  /////////////////////////////////////
std::vector<std::vector<int32_t>> GenerateRandomMatrixInt32(size_t n, int32_t maxValue)
{
    OpenSSL_URBG rng;
    std::uniform_int_distribution<int32_t> dist(0, maxValue);

    std::vector<std::vector<int32_t>> matrix(n, std::vector<int32_t>(n));

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            matrix[i][j] = dist(rng);

    return matrix;
}

/////////////////////////////////////   Secret s  /////////////////////////////////////
int32_t sample_eta_centered_binomial(uint8_t eta)
{
    OpenSSL_URBG rng;
    std::uniform_int_distribution<uint8_t> dis(0, 1);
    uint8_t sum1 = 0, sum2 = 0;
    for (uint8_t i = 0; i < eta; ++i)
    {
        sum1 += dis(rng);
        sum2 += dis(rng);
    }
    return (int32_t)sum1 - (int32_t)sum2;
}

std::vector<int32_t> sample_vector_binomial(uint32_t n)
{
    std::vector<int32_t> result;
    result.reserve(n);
    for (uint32_t i = 0; i < n; ++i)
        result.push_back(sample_eta_centered_binomial(3));
    return result;
}

/////////////////////////////////////   Gaussian noise  /////////////////////////////////////
int32_t sample_discrete_gaussian(double sigma)
{
    OpenSSL_URBG rng;
    std::normal_distribution<double> norm(0.0, sigma);

    double x = norm(rng);

    double bound = 6.0 * sigma;
    if (x > bound)  x = bound;
    if (x < -bound) x = -bound;

    long long k = std::llround(x);

    if (k > std::numeric_limits<int32_t>::max()) k = std::numeric_limits<int32_t>::max();
    if (k < std::numeric_limits<int32_t>::min()) k = std::numeric_limits<int32_t>::min();
    return static_cast<int32_t>(k);
}

std::vector<int32_t> GenerateGaussianVector(size_t n)
{
    std::vector<int32_t> vec;
    vec.reserve(n);
    for (size_t i = 0; i < n; ++i)
        vec.push_back(sample_discrete_gaussian(2.3));
    return vec;
}

int32_t getRandomInt(int min, int max)
{
    OpenSSL_URBG rng;
    std::uniform_int_distribution<> distr(min, max);
    return distr(rng);
}

int random_bit()
{
    OpenSSL_URBG urbg;
    std::uniform_int_distribution<int> bit(0, 1);
    return bit(urbg);
}

/////////////////////////////////////   Helpers  /////////////////////////////////////
std::vector<int32_t> concat(const std::vector<int32_t> &a, const std::vector<int32_t> &b)
{
    std::vector<int32_t> out;
    out.reserve(a.size() + b.size());
    for (size_t i = 0; i < a.size(); ++i)
        out.push_back(a[i]);
    for (size_t i = 0; i < b.size(); ++i)
        out.push_back(b[i]);
    return out;
}

std::vector<int32_t> flatten_matrix(const std::vector<std::vector<int32_t>> &A)
{
    std::vector<int32_t> out;
    if (A.empty()) return out;
    const size_t n = A.size();
    out.reserve(n * A[0].size());
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            out.push_back(A[i][j]);
    return out;
}

void transposeA(const std::vector<std::vector<int32_t>> &A, std::vector<std::vector<int32_t>> &AT)
{
    const size_t n = A.size();
    AT.assign(n, std::vector<int32_t>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            AT[j][i] = A[i][j];
}