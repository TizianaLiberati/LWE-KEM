#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

// prende 1 byte (LSB) per ogni int32_t
std::vector<uint8_t> int32_to_bytes_keccak(const std::vector<int32_t>& v);

// SHAKE256 one-shot: input bytes -> output bytes (lunghezza scelta)
std::vector<uint8_t> shake256_keccak(const std::vector<uint8_t>& input, size_t out_len_bytes);

// seed(int32_t) -> bytes -> SHAKE256 per ottenere coins
std::vector<uint8_t> xof_coins_keccak(const std::vector<int32_t>& seed_int32,
                                     size_t n,
                                     size_t msg_bits);

// SHA3-256: input vector<int32_t> -> digest (32 byte) come vector<int32_t> [0..255]
std::vector<int32_t> SHA3_256_keccak(const std::vector<int32_t>& in);
