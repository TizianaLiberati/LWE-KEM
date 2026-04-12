#pragma once
#include <cstddef>
#include <vector>
#include <cstdint>

// ============================================================================
// CORE HASH FUNCTIONS
// ============================================================================

std::vector<uint8_t> int32_to_bytes_openssl(const std::vector<int32_t>& v);

std::vector<uint8_t> shake256_openssl(const std::vector<uint8_t>& input, 
                                     size_t out_len_bytes);

std::vector<uint8_t> xof_coins_openssl(const std::vector<int32_t>& seed_int32, 
                                        size_t n, 
                                        size_t msg_bits);

std::vector<int32_t> SHA3_256_openssl(const std::vector<int32_t>& in);

// ============================================================================
// OPTIMIZED INTERFACES
// ============================================================================

// In-place byte conversion to avoid allocation
void int32_to_bytes_inplace(const std::vector<int32_t>& v, 
                            uint8_t* out, 
                            size_t out_len);

// Streaming XOF with callback for memory-efficient processing
void shake256_openssl_streaming(const std::vector<uint8_t>& input,
                                 size_t total_out_bytes,
                                 size_t chunk_size,
                                 void (*callback)(const uint8_t*, size_t, void*),
                                 void* user_data);

// Buffered XOF that reuses output buffer
void xof_coins_openssl_buffered(const std::vector<int32_t>& seed_int32,
                                  size_t n,
                                  size_t msg_bits,
                                  std::vector<uint8_t>& out_buffer);

// In-place SHA3-256 for raw byte arrays
void SHA3_256_openssl_inplace(const uint8_t* in, 
                               size_t in_len, 
                               uint8_t out[32]);

