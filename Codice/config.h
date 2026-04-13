#pragma once
// config.h
// Central definition of all runtime-tunable LWE-KEM parameters.
//
// Design principles:
//   - Every parameter that affects correctness OR performance lives here.
//   - Defaults are the mathematically safe values for n up to 16384, q=3329.
//   - The struct is plain-old-data; no OpenMP state is stored here.
//   - Thread safety: the Config object is written once (during parsing) before
//     any OpenMP region is entered, then read-only for the rest of the run.
// ============================================================================

#include <cstdint>
#include <string>

// ---------------------------------------------------------------------------
// Noise distribution selector
// ---------------------------------------------------------------------------
// BINOMIAL_ETA1 : X = b1 - b2, b_i ~ Bernoulli(0.5), range {-1,0,1}, σ=0.707
//                 Required for n > ~400 with q=3329 (see noise analysis).
// BINOMIAL_ETA3 : X = Σ3(b_i) - Σ3(b_i), range {-3,...,3}, σ=1.22
//                 Safe only for n ≤ 256 with q=3329.
// GAUSSIAN      : X ~ N(0,σ²), rounded to integer.
//                 NOT safe for any n with q=3329 (σ=2.3 → 6σ·√n >> q/4).
//                 Kept for research/comparison only.
enum class NoiseKind { BINOMIAL_ETA1, BINOMIAL_ETA3, GAUSSIAN, UNIFORM_CENTERED_3 };

inline const char* noise_kind_str(NoiseKind k) {
    switch (k) {
        case NoiseKind::BINOMIAL_ETA1: return "BINOMIAL_ETA1";
        case NoiseKind::BINOMIAL_ETA3: return "BINOMIAL_ETA3";
        case NoiseKind::GAUSSIAN:      return "GAUSSIAN";
        case NoiseKind::UNIFORM_CENTERED_3: return "UNIFORM_CENTERED_3";
    }
    return "UNKNOWN";
}

// ---------------------------------------------------------------------------
// LweKemConfig — all parameters in one struct
// ---------------------------------------------------------------------------
struct LweKemConfig {

    // ---- Problem parameters ------------------------------------------------
    uint32_t q  = 3329;    // modulus  (Kyber default; must be prime)
    uint32_t n  = 256;     // lattice dimension
    int      N  = 10;      // number of benchmark iterations

    // ---- Secret / error distributions -------------------------------------
    // s_eta : eta for the secret key s (Binomial).  3 = Kyber-standard.
    // e_kind, r_kind, e1_kind, e2_kind: distribution for each noise term.
    //
    // CORRECTNESS CONSTRAINT (automatically checked at startup):
    //   total_noise_6sigma < q/4
    //   where total_noise_6sigma = 6 * sqrt(n*(σ_e²·σ_r²  +  σ_s²·σ_e1²) + σ_e2²)
    //   With q=3329, n=16384: only BINOMIAL_ETA1 for e/r/e1/e2 is safe.
    uint32_t   s_eta    = 3;                   // secret key eta
    NoiseKind  e_kind   = NoiseKind::BINOMIAL_ETA1; // keygen error e
    NoiseKind  r_kind   = NoiseKind::BINOMIAL_ETA1; // encaps randomness r
    NoiseKind  e1_kind  = NoiseKind::BINOMIAL_ETA1; // encaps error e1
    NoiseKind  e2_kind  = NoiseKind::BINOMIAL_ETA1; // encaps scalar e2
    double     gauss_sigma = 2.3;              // sigma if NoiseKind==GAUSSIAN

    // ---- Optional reproducibility seed ------------------------------------
    // If seed == 0  → use system entropy (non-reproducible, default).
    // If seed != 0  → passed to OpenSSL RAND_seed before any RNG use,
    //                 making the run deterministic for debugging.
    uint64_t seed = 0;

    // ---- OpenMP tuning ----------------------------------------------------
    // num_threads == 0 means "use OMP_NUM_THREADS / hardware default".
    int num_threads = 0;

    // ---- Output control ---------------------------------------------------
    bool verbose       = true;   // print parameter table before CSV
    bool warn_mismatch = true;   // print WARNING to stderr if mismatches > 0
    bool suggest_fix   = true;   // suggest parameter changes when mismatches > 0
};
