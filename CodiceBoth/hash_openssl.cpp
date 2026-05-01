#include "hash_openssl.h"

#include <openssl/evp.h>
#include <openssl/err.h>

#include <stdexcept>
#include <string>

// ---- helper errore OpenSSL
static void throw_ossl(const char* where)
{
    unsigned long e = ERR_get_error();
    char buf[256];
    ERR_error_string_n(e, buf, sizeof(buf));
    throw std::runtime_error(std::string(where) + ": " + buf);
}

std::vector<uint8_t> int32_to_bytes_openssl(const std::vector<int32_t>& v)
{
    // prende il byte basso
    std::vector<uint8_t> out;
    out.reserve(v.size());
    for (int32_t x : v)
        out.push_back(static_cast<uint8_t>(x & 0xFF));
    return out;
}

std::vector<uint8_t> shake256_openssl(const std::vector<uint8_t>& input, size_t out_len_bytes)
{
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) throw_ossl("EVP_MD_CTX_new");

    const EVP_MD* md = EVP_shake256();
    if (!md) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("EVP_shake256() returned null (SHAKE256 not available?)");
    }

    if (EVP_DigestInit_ex(ctx, md, nullptr) != 1) { EVP_MD_CTX_free(ctx); throw_ossl("EVP_DigestInit_ex(SHAKE256)"); }
    if (EVP_DigestUpdate(ctx, input.data(), input.size()) != 1) { EVP_MD_CTX_free(ctx); throw_ossl("EVP_DigestUpdate(SHAKE256)"); }

    std::vector<uint8_t> out(out_len_bytes);

    // XOF output: EVP_DigestFinalXOF
    if (EVP_DigestFinalXOF(ctx, out.data(), out.size()) != 1) { EVP_MD_CTX_free(ctx); throw_ossl("EVP_DigestFinalXOF(SHAKE256)"); }

    EVP_MD_CTX_free(ctx);
    return out;
    // EVP_DigestFinalXOF è l’API corretta per SHAKE XOF. :contentReference[oaicite:1]{index=1}
}

std::vector<uint8_t> xof_coins_openssl(const std::vector<int32_t>& seed_int32, size_t n, size_t msg_bits)
{
    std::vector<uint8_t> seed_bytes = int32_to_bytes_openssl(seed_int32);

    size_t total_bytes_per_bit = 2 * n + 1;
    size_t out_bytes = msg_bits * total_bytes_per_bit;

    return shake256_openssl(seed_bytes, out_bytes);
}

std::vector<int32_t> SHA3_256_openssl(const std::vector<int32_t>& in)
{
    // conversione: ogni int32 -> 4 byte little-endian
    std::vector<uint8_t> bytes;
    bytes.reserve(in.size() * 4);

    for (int32_t x : in) {
        uint32_t ux = static_cast<uint32_t>(x);
        bytes.push_back(static_cast<uint8_t>(ux & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 8)  & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 16) & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 24) & 0xFF));
    }

    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) throw_ossl("EVP_MD_CTX_new");

    const EVP_MD* md = EVP_sha3_256();
    if (!md) { EVP_MD_CTX_free(ctx); throw std::runtime_error("EVP_sha3_256() returned null (SHA3-256 not available?)"); }

    if (EVP_DigestInit_ex(ctx, md, nullptr) != 1) { EVP_MD_CTX_free(ctx); throw_ossl("EVP_DigestInit_ex(SHA3-256)"); }
    if (EVP_DigestUpdate(ctx, bytes.data(), bytes.size()) != 1) { EVP_MD_CTX_free(ctx); throw_ossl("EVP_DigestUpdate(SHA3-256)"); }

    unsigned char digest[32];
    unsigned int digest_len = 0;

    if (EVP_DigestFinal_ex(ctx, digest, &digest_len) != 1) { EVP_MD_CTX_free(ctx); throw_ossl("EVP_DigestFinal_ex(SHA3-256)"); }
    EVP_MD_CTX_free(ctx);

    if (digest_len != 32) throw std::runtime_error("SHA3-256 digest length != 32");

    std::vector<int32_t> out(32);
    for (size_t i = 0; i < 32; ++i)
        out[i] = static_cast<int32_t>(digest[i]);

    return out;
    // EVP_* è l’interfaccia raccomandata per digest in OpenSSL. :contentReference[oaicite:2]{index=2}
}