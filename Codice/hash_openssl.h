#pragma once
#include <vector>
#include <cstdint>

std::vector<uint8_t> int32_to_bytes_openssl(const std::vector<int32_t>& v);

std::vector<uint8_t> shake256_openssl(const std::vector<uint8_t>& input, size_t out_len_bytes);

std::vector<uint8_t> xof_coins_openssl(const std::vector<int32_t>& seed_int32, size_t n, size_t msg_bits);

std::vector<int32_t> SHA3_256_openssl(const std::vector<int32_t>& in);