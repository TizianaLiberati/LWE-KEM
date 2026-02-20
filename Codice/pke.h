#pragma once
#include <vector>
#include <cstdint>
#include <omp.h>

#include "utils.h"

void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &s_k, std::vector<int32_t> &t);

void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i, 
    std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2, const std::vector<std::vector<int32_t>> &AT);

void Decrypt(int32_t v_i, const std::vector<int32_t> &u, const std::vector<int32_t> &s_k, uint32_t q, int32_t &decrypt_i);