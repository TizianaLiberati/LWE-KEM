#pragma once

#include <openssl/rand.h>
#include <openssl/err.h>

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

// #pragma acc routine seq
inline void openssl_randomize(std::uint8_t* out, size_t len)
{
    if (RAND_bytes(out, static_cast<int>(len)) != 1) {
        unsigned long e = ERR_get_error();
        char buf[256];
        ERR_error_string_n(e, buf, sizeof(buf));
        throw std::runtime_error(std::string("RAND_bytes failed: ") + buf);
    }
}

// Adattatore URBG per <random>
struct OpenSSL_URBG
{
    using result_type = std::uint32_t;

    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    result_type operator()()
    {
        result_type x;
        openssl_randomize(reinterpret_cast<std::uint8_t*>(&x), sizeof(x));
        return x;
    }
};