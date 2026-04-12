#include "hash_openssl.h"

#include <openssl/evp.h>
#include <openssl/err.h>

#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>

// ============================================================================
// ERROR HANDLING
// ============================================================================

static void throw_ossl(const char* where) {
    unsigned long e = ERR_get_error();
    char buf[256];
    ERR_error_string_n(e, buf, sizeof(buf));
    throw std::runtime_error(std::string(where) + ": " + buf);
}

// ============================================================================
// OPTIMIZED BYTE CONVERSION
// ============================================================================

std::vector<uint8_t> int32_to_bytes_openssl(const std::vector<int32_t>& v) {
    // OPTIMIZED: Pre-allocate and use pointer arithmetic
    std::vector<uint8_t> out;
    out.reserve(v.size());
    
    const int32_t* data = v.data();
    const size_t len = v.size();
    
    for (size_t i = 0; i < len; ++i) {
        out.push_back(static_cast<uint8_t>(data[i] & 0xFF));
    }
    
    return out;
}

// OPTIMIZED: Direct conversion with fixed-size output
void int32_to_bytes_inplace(const std::vector<int32_t>& v, 
                            uint8_t* out, 
                            size_t out_len) {
    const size_t copy_len = std::min(v.size(), out_len);
    const int32_t* data = v.data();
    
    for (size_t i = 0; i < copy_len; ++i) {
        out[i] = static_cast<uint8_t>(data[i] & 0xFF);
    }
}

// ============================================================================
// STREAMING SHAKE256 XOF
// ============================================================================

// SHAKE256 context for incremental output
class Shake256Ctx {
public:
    EVP_MD_CTX* ctx;
    
    Shake256Ctx() {
        ctx = EVP_MD_CTX_new();
        if (!ctx) throw_ossl("EVP_MD_CTX_new");
        
        const EVP_MD* md = EVP_shake256();
        if (!md) {
            EVP_MD_CTX_free(ctx);
            throw std::runtime_error("EVP_shake256() returned null");
        }
        
        if (EVP_DigestInit_ex(ctx, md, nullptr) != 1) {
            throw_ossl("EVP_DigestInit_ex(SHAKE256)");
        }
    }
    
    ~Shake256Ctx() {
        if (ctx) EVP_MD_CTX_free(ctx);
    }
    
    void update(const uint8_t* data, size_t len) {
        if (EVP_DigestUpdate(ctx, data, len) != 1) {
            throw_ossl("EVP_DigestUpdate(SHAKE256)");
        }
    }
    
    void squeeze(uint8_t* out, size_t len) {
        // XOF output: EVP_DigestFinalXOF
        if (EVP_DigestFinalXOF(ctx, out, len) != 1) {
            throw_ossl("EVP_DigestFinalXOF(SHAKE256)");
        }
    }
    
    // For additional output after finalization, reinitialize
    void reset() {
        const EVP_MD* md = EVP_shake256();
        if (EVP_DigestInit_ex(ctx, md, nullptr) != 1) {
            throw_ossl("EVP_DigestInit_ex(SHAKE256 reset)");
        }
    }
    
    // Disable copy
    Shake256Ctx(const Shake256Ctx&) = delete;
    Shake256Ctx& operator=(const Shake256Ctx&) = delete;
};

// ============================================================================
// OPTIMIZED SHAKE256 OUTPUT
// ============================================================================

std::vector<uint8_t> shake256_openssl(const std::vector<uint8_t>& input, 
                                     size_t out_len_bytes) {
    Shake256Ctx shake;
    shake.update(input.data(), input.size());
    
    std::vector<uint8_t> out(out_len_bytes);
    shake.squeeze(out.data(), out_len_bytes);
    
    return out;
}

// OPTIMIZED: Streaming XOF that processes output in chunks
// Reduces memory allocation for large outputs
void shake256_openssl_streaming(const std::vector<uint8_t>& input,
                                 size_t total_out_bytes,
                                 size_t chunk_size,
                                 void (*callback)(const uint8_t*, size_t, void*),
                                 void* user_data) {
    Shake256Ctx shake;
    shake.update(input.data(), input.size());
    
    std::vector<uint8_t> chunk(chunk_size);
    size_t remaining = total_out_bytes;
    
    while (remaining > 0) {
        size_t to_read = std::min(chunk_size, remaining);
        shake.squeeze(chunk.data(), to_read);
        callback(chunk.data(), to_read, user_data);
        remaining -= to_read;
    }
}

// ============================================================================
// OPTIMIZED XOF COINS GENERATION
// ============================================================================

std::vector<uint8_t> xof_coins_openssl(const std::vector<int32_t>& seed_int32, 
                                        size_t n, 
                                        size_t msg_bits) {
    // OPTIMIZED: Stack allocation for small seeds, avoid heap allocation
    uint8_t seed_bytes[64];  // Assume max 64 bytes for seed
    size_t seed_len = std::min(seed_int32.size(), size_t(64));
    int32_to_bytes_inplace(seed_int32, seed_bytes, seed_len);
    
    size_t total_bytes_per_bit = 2 * n + 1;
    size_t out_bytes = msg_bits * total_bytes_per_bit;
    
    // OPTIMIZED: Direct vector construction
    std::vector<uint8_t> result;
    result.reserve(out_bytes);
    
    Shake256Ctx shake;
    shake.update(seed_bytes, seed_len);
    
    // Allocate result and squeeze directly
    result.resize(out_bytes);
    shake.squeeze(result.data(), out_bytes);
    
    return result;
}

// OPTIMIZED: Thread-safe version with pre-allocated buffer
void xof_coins_openssl_buffered(const std::vector<int32_t>& seed_int32,
                                  size_t n,
                                  size_t msg_bits,
                                  std::vector<uint8_t>& out_buffer) {
    uint8_t seed_bytes[64];
    size_t seed_len = std::min(seed_int32.size(), size_t(64));
    int32_to_bytes_inplace(seed_int32, seed_bytes, seed_len);
    
    size_t total_bytes_per_bit = 2 * n + 1;
    size_t out_bytes = msg_bits * total_bytes_per_bit;
    
    // Resize output buffer (reuses memory if already large enough)
    out_buffer.resize(out_bytes);
    
    Shake256Ctx shake;
    shake.update(seed_bytes, seed_len);
    shake.squeeze(out_buffer.data(), out_bytes);
}

// ============================================================================
// OPTIMIZED SHA3-256
// ============================================================================

std::vector<int32_t> SHA3_256_openssl(const std::vector<int32_t>& in) {
    // OPTIMIZED: Stack buffer for conversion, avoid heap allocation for small inputs
    const size_t max_stack_size = 256;
    uint8_t stack_buffer[max_stack_size * 4];
    
    size_t in_len = in.size();
    size_t bytes_len = in_len * 4;
    
    uint8_t* bytes = (bytes_len <= max_stack_size * 4) ? stack_buffer 
                                                          : new uint8_t[bytes_len];
    
    // Fast little-endian conversion
    const int32_t* in_data = in.data();
    for (size_t i = 0; i < in_len; ++i) {
        uint32_t ux = static_cast<uint32_t>(in_data[i]);
        size_t base = i * 4;
        bytes[base]     = static_cast<uint8_t>(ux & 0xFF);
        bytes[base + 1] = static_cast<uint8_t>((ux >> 8)  & 0xFF);
        bytes[base + 2] = static_cast<uint8_t>((ux >> 16) & 0xFF);
        bytes[base + 3] = static_cast<uint8_t>((ux >> 24) & 0xFF);
    }
    
    // SHA3-256 computation
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) {
        if (bytes != stack_buffer) delete[] bytes;
        throw_ossl("EVP_MD_CTX_new");
    }
    
    const EVP_MD* md = EVP_sha3_256();
    if (!md) {
        EVP_MD_CTX_free(ctx);
        if (bytes != stack_buffer) delete[] bytes;
        throw std::runtime_error("EVP_sha3_256() returned null");
    }
    
    if (EVP_DigestInit_ex(ctx, md, nullptr) != 1) {
        EVP_MD_CTX_free(ctx);
        if (bytes != stack_buffer) delete[] bytes;
        throw_ossl("EVP_DigestInit_ex(SHA3-256)");
    }
    
    if (EVP_DigestUpdate(ctx, bytes, bytes_len) != 1) {
        EVP_MD_CTX_free(ctx);
        if (bytes != stack_buffer) delete[] bytes;
        throw_ossl("EVP_DigestUpdate(SHA3-256)");
    }
    
    unsigned char digest[32];
    unsigned int digest_len = 0;
    
    if (EVP_DigestFinal_ex(ctx, digest, &digest_len) != 1) {
        EVP_MD_CTX_free(ctx);
        if (bytes != stack_buffer) delete[] bytes;
        throw_ossl("EVP_DigestFinal_ex(SHA3-256)");
    }
    
    EVP_MD_CTX_free(ctx);
    if (bytes != stack_buffer) delete[] bytes;
    
    if (digest_len != 32) {
        throw std::runtime_error("SHA3-256 digest length != 32");
    }
    
    // OPTIMIZED: Direct vector initialization
    std::vector<int32_t> out(32);
    for (size_t i = 0; i < 32; ++i) {
        out[i] = static_cast<int32_t>(digest[i]);
    }
    
    return out;
}

// OPTIMIZED: In-place SHA3-256 for known-size inputs
void SHA3_256_openssl_inplace(const uint8_t* in, 
                               size_t in_len, 
                               uint8_t out[32]) {
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) throw_ossl("EVP_MD_CTX_new");
    
    const EVP_MD* md = EVP_sha3_256();
    if (!md) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("EVP_sha3_256() returned null");
    }
    
    if (EVP_DigestInit_ex(ctx, md, nullptr) != 1) {
        EVP_MD_CTX_free(ctx);
        throw_ossl("EVP_DigestInit_ex(SHA3-256)");
    }
    
    if (EVP_DigestUpdate(ctx, in, in_len) != 1) {
        EVP_MD_CTX_free(ctx);
        throw_ossl("EVP_DigestUpdate(SHA3-256)");
    }
    
    unsigned int digest_len = 0;
    if (EVP_DigestFinal_ex(ctx, out, &digest_len) != 1) {
        EVP_MD_CTX_free(ctx);
        throw_ossl("EVP_DigestFinal_ex(SHA3-256)");
    }
    
    EVP_MD_CTX_free(ctx);
}

