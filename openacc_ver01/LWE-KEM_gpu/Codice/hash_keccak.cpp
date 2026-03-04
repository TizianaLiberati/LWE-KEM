#include "hash_keccak.h"

#include <stdexcept>
#include <string>

// KeccakCodePackage/XKCP è C → evita name mangling
extern "C" {
#include "SimpleFIPS202.h"  // for SHA3_256(), SHAKE256()
}

// 1 byte per int32 (LSB)
std::vector<uint8_t> int32_to_bytes_keccak(const std::vector<int32_t>& v)
{
    std::vector<uint8_t> out;
    out.reserve(v.size());
    for (int32_t x : v) {
        out.push_back(static_cast<uint8_t>(x & 0xFF));
    }
    return out;
}

std::vector<uint8_t> shake256_keccak(const std::vector<uint8_t>& input, size_t out_len_bytes)
{
    std::vector<uint8_t> out(out_len_bytes);

    // API SimpleFIPS202: int SHAKE256(out, outLen, in, inLen)
    const int rc = SHAKE256(out.data(), out.size(), input.data(), input.size());
    if (rc != 0) {
        throw std::runtime_error("SHAKE256() failed (Keccak/XKCP)");
    }
    return out;
}

std::vector<uint8_t> xof_coins_keccak(const std::vector<int32_t>& seed_int32, size_t n, size_t msg_bits)
{
    // seed -> bytes
    std::vector<uint8_t> seed_bytes = int32_to_bytes_keccak(seed_int32);

    // per ogni bit: (2*n + 1) byte
    const size_t total_bytes_per_bit = 2 * n + 1;
    const size_t out_bytes = msg_bits * total_bytes_per_bit;

    return shake256_keccak(seed_bytes, out_bytes);
}

std::vector<int32_t> SHA3_256_keccak(const std::vector<int32_t>& in)
{
    // ogni int32_t viene serializzato in 4 byte LE
    std::vector<uint8_t> bytes;
    bytes.reserve(in.size() * 4);

    for (int32_t x : in) {
        uint32_t ux = static_cast<uint32_t>(x);
        bytes.push_back(static_cast<uint8_t>(ux & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 8)  & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 16) & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 24) & 0xFF));
    }

    std::vector<uint8_t> digest(32);
    const int rc = SHA3_256(digest.data(), bytes.data(), bytes.size());
    if (rc != 0) {
        throw std::runtime_error("SHA3_256() failed (Keccak/XKCP)");
    }

    // output: 32 elementi [0..255] in int32_t
    std::vector<int32_t> out(digest.size());
    for (size_t i = 0; i < digest.size(); ++i)
        out[i] = static_cast<int32_t>(digest[i]);

    return out;
}