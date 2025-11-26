#pragma once
#include <vector>
#include <cstdint>

// Creiamo una struttura per i noise che useremo spesso
struct NoiseTriple
{
    std::vector<int32_t> r;
    std::vector<int32_t> e1;
    int32_t e2;
};

int32_t sample_e2(uint8_t coins);

double uniform01_from_byte(uint8_t b);

double normal_cdf(double x);

double normal_icdf(double u);

double sample_normal01_from_byte(uint8_t coins);

int32_t sample_discrete_gaussian_sigma_from_byte(uint8_t coins, double sigma);

NoiseTriple GenerateNoisesForOneBit(const std::vector<uint8_t>& coins, size_t& pos, const size_t n);