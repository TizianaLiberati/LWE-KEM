#pragma once
#include <vector>
#include <cstdint>
#include <omp.h>

#include "utils.h"
#include "hash_openssl.h"
#include "noise.h"
#include "pke.h"
#include "rng_openssl.h"

void Encaps(uint32_t n, uint32_t q, 
            std::vector<int32_t> &t, 
            std::vector<int32_t> &c, 
            const std::vector<std::vector<int32_t>> &A, 
            const std::vector<std::vector<int32_t>> &AT, 
            std::vector<int32_t> &Hash_K);

void Decaps(uint32_t n, uint32_t q, 
            const std::vector<int32_t> &t, 
            const std::vector<int32_t> &s_k, 
            const std::vector<int32_t> &c, 
            std::vector<int32_t> &Hash_K, 
            const std::vector<std::vector<int32_t>> &A, 
            const std::vector<std::vector<int32_t>> &AT);

