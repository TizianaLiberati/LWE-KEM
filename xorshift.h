/*  xorshift.h  –  device-side PRNG for OpenMP target offloading kernels
 *
 *  NOTE: These functions are NOT wrapped in #pragma omp declare target
 *  here because nvc++ 24.5 hits an ICE when declare-target inline
 *  functions are included from a header.  Instead, pke.cpp embeds
 *  copies of these functions inside its own declare-target block.
 *
 *  This header exists for CPU-only callers (e.g. kem.cpp).
 */
#pragma once
#include <cstdint>
#include <cmath>

/* ------------------------------------------------------------------ */
/*  Core xorshift64*                                                   */
/* ------------------------------------------------------------------ */
static inline uint64_t xorshift64star(uint64_t *state)
{
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

/* Uniform integer in [0, max_val] */
static inline int32_t xorshift_uniform(uint64_t *state, int32_t max_val)
{
    return (int32_t)(xorshift64star(state) % (uint64_t)(max_val + 1));
}

/* ------------------------------------------------------------------ */
/*  Sampling helpers                                                   */
/* ------------------------------------------------------------------ */

/* Centered-binomial  η = eta  →  value in [-eta, +eta] */
static inline int32_t device_sample_cbd(uint64_t *state, int eta)
{
    int32_t a = 0, b = 0;
    for (int i = 0; i < eta; ++i) {
        a += xorshift_uniform(state, 1);
        b += xorshift_uniform(state, 1);
    }
    return a - b;
}

/* Discrete Gaussian  N(0, σ²)  via Box-Muller, rounded to int32      */
static inline int32_t device_sample_gaussian(uint64_t *state, double sigma)
{
    double u1 = ((double)(xorshift64star(state) % 1000000000ULL) + 1.0) / 1000000001.0;
    double u2 = ((double)(xorshift64star(state) % 1000000000ULL) + 1.0) / 1000000001.0;

    double z = sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
    double val = sigma * z;

    double bnd = 6.0 * sigma;
    if (val >  bnd) val =  bnd;
    if (val < -bnd) val = -bnd;

    return (int32_t)(val >= 0.0 ? val + 0.5 : val - 0.5);
}
