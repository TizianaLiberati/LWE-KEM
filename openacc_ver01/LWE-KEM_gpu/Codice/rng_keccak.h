#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

#if defined(_WIN32)
  #include <windows.h>
  #include <bcrypt.h>
  #pragma comment(lib, "bcrypt.lib")
#elif defined(__APPLE__)
  #include <stdlib.h> // arc4random_buf
#else
  #include <errno.h>
  #include <sys/random.h> // getrandom (Linux)
  #include <unistd.h>
#endif

extern "C" {
#include "SimpleFIPS202.h" // SHAKE256()
}

// --------- OS entropy helper ----------
inline void keccak_os_random(uint8_t* out, size_t len)
{
#if defined(_WIN32)
    if (BCryptGenRandom(nullptr, out, (ULONG)len, BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
        throw std::runtime_error("BCryptGenRandom failed");
    }
#elif defined(__APPLE__)
    arc4random_buf(out, len);
#else
    // Linux: getrandom()
    size_t off = 0;
    while (off < len) {
        ssize_t r = getrandom(out + off, len - off, 0);
        if (r < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error("getrandom failed");
        }
        off += (size_t)r;
    }
#endif
}

// --------- Simple SHAKE256-DRBG (seeded once) ----------
class KeccakShake256DRBG {
public:
    KeccakShake256DRBG()
    {
        keccak_os_random(key_.data(), key_.size());
        counter_ = 0;
        used_ = buffer_.size(); // forza refill al primo uso
    }

    void randomize(uint8_t* out, size_t len)
    {
        std::lock_guard<std::mutex> lock(mu_);
        while (len > 0) {
            if (used_ == buffer_.size()) refill_();

            const size_t avail = buffer_.size() - used_;
            const size_t take = (len < avail) ? len : avail;

            std::memcpy(out, buffer_.data() + used_, take);
            used_ += take;

            out += take;
            len -= take;
        }
    }

private:
    void refill_()
    {
        // input = key || counter
        std::array<uint8_t, 32 + 8> in{};
        std::memcpy(in.data(), key_.data(), 32);
        for (int i = 0; i < 8; ++i) {
            in[32 + i] = (uint8_t)((counter_ >> (8 * i)) & 0xFF);
        }
        counter_++;

        const int rc = SHAKE256(buffer_.data(), buffer_.size(), in.data(), in.size());
        if (rc != 0) throw std::runtime_error("SHAKE256 failed in DRBG");
        used_ = 0;
    }

    std::array<uint8_t, 32> key_{};
    uint64_t counter_{0};

    std::array<uint8_t, 256> buffer_{}; // buffer interno
    size_t used_{0};

    std::mutex mu_;
};

inline KeccakShake256DRBG& global_rng_keccak()
{
    static KeccakShake256DRBG rng;
    return rng;
}

// URBG adapter “stile Botan_URBG”
struct Keccak_URBG
{
    using result_type = std::uint32_t;

    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    result_type operator()()
    {
        result_type x = 0;
        global_rng_keccak().randomize(reinterpret_cast<std::uint8_t*>(&x), sizeof(x));
        return x;
    }
};
