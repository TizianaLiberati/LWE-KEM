// Codice_modificato.cpp
// LWE-KEM benchmark driver — config-file driven.
//
// Usage:
//   ./lwe_kem config.txt          <- preferred: all parameters from file
//   ./lwe_kem N n                 <- legacy CLI fallback (N iterations, dimension n)
//
// OpenMP notes:
//   omp_set_max_active_levels(1) prevents nested thread-team explosion.
//   Parallelism lives inside KeyGen, Encaps, Decaps only.
//   The outer N-iteration loop is intentionally sequential (per-iteration timing).
// ============================================================================

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>
#include <openssl/rand.h>

#include "config.h"
#include "config_parser.h"
#include "kem.h"
#include "noise.h"
#include "pke.h"
#include "utils.h"
#include "hash_openssl.h"
#include "utils.h"

// ============================================================================
// apply_seed
// ============================================================================
// Mixes the user-supplied seed into OpenSSL's PRNG before any RNG use.
// Provides approximate reproducibility for debugging; a truly deterministic
// run would require a custom DRBG, which is outside this KEM's scope.

static void apply_seed(uint64_t seed) {
    if (seed == 0) return;
    unsigned char buf[8];
    for (int i = 0; i < 8; ++i)
        buf[i] = static_cast<unsigned char>((seed >> (8*i)) & 0xFF);
    RAND_seed(buf, 8);
    std::cerr << "[seed] OpenSSL PRNG seeded with " << seed
              << " (approximate reproducibility)\n";
}

// ============================================================================
// looks_like_file
// ============================================================================
// Heuristic: if the single argument is not purely numeric, treat it as a
// filename so that "./lwe_kem 10 256" keeps working as legacy CLI.

static bool looks_like_file(const char *arg) {
    for (const char *p = arg; *p; ++p)
        if (!std::isdigit(static_cast<unsigned char>(*p)))
            return true;
    return false;
}

// ============================================================================
// main
// ============================================================================

int main(int argc, char **argv) {

    // Disable nested OpenMP teams before any parallel work.
    omp_set_max_active_levels(1);

    // ---- Load configuration ------------------------------------------------
    LweKemConfig cfg;
    std::string  config_source = "built-in defaults";

    try {
        if (argc == 2 && looks_like_file(argv[1])) {
            // ./lwe_kem config.txt
            cfg = parse_config(argv[1]);
            config_source = argv[1];

        } else if (argc >= 3 && !looks_like_file(argv[1])) {
            // ./lwe_kem N n   (legacy)
            cfg = parse_cli_fallback(argc, argv);
            config_source = "CLI arguments";

        } else if (argc == 1) {
            // No arguments: try default config file in CWD.
            std::ifstream probe("config.txt");
            if (probe.good()) {
                probe.close();
                cfg = parse_config("config.txt");
                config_source = "config.txt (auto-detected)";
                std::cerr << "[config] No argument — loaded config.txt\n";
            } else {
                std::cerr << "[config] No argument and no config.txt found; "
                             "using built-in defaults (N=10, n=256).\n";
            }
        } else {
            std::cerr << "Usage:\n"
                      << "  " << argv[0] << " config.txt\n"
                      << "  " << argv[0] << " N n\n";
            return 1;
        }
    } catch (const std::exception &ex) {
        std::cerr << "Configuration error: " << ex.what() << "\n";
        return 1;
    }

    // ---- Apply thread count ------------------------------------------------
    if (cfg.num_threads > 0)
        omp_set_num_threads(cfg.num_threads);

    // ---- Apply seed --------------------------------------------------------
    apply_seed(cfg.seed);

    // ---- Print configuration -----------------------------------------------
    if (cfg.verbose) {
        std::cout << "[config source: " << config_source << "]\n\n";
        print_config(cfg, std::cout);
        std::cout << "\n";
    }

    // ---- Validate noise parameters -----------------------------------------
    bool noise_ok = validate_config(cfg, std::cout);
    std::cout << "\n";

    if (!noise_ok && cfg.warn_mismatch) {
        std::cerr << "WARNING: noise parameters are incompatible with "
                     "n=" << cfg.n << ", q=" << cfg.q << ".\n"
                  << "         Decryption mismatches are EXPECTED on this run.\n"
                  << "         See suggestions above, or add to config.txt:\n"
                  << "           noise_kind = binomial_eta1\n\n";
    }

    // ---- Extract scalar parameters -----------------------------------------
    const uint32_t q = cfg.q;
    const uint32_t n = cfg.n;
    const int      N = cfg.N;

    // Print run header.
    std::cout << "LWE-KEM Benchmark\n"
              << "=================\n"
              << "N=" << N
              << "  n=" << n
              << "  q=" << q
              << "  threads=" << omp_get_max_threads()
              << "\n\n";

    // ---- Benchmark loop ----------------------------------------------------
    long long sum_keygen_us = 0;
    long long sum_encaps_us = 0;
    long long sum_decaps_us = 0;
    int       mismatches    = 0;

    const auto wall_start = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k) {

        // ---- KeyGen --------------------------------------------------------
        std::vector<std::vector<int32_t>> A;
        std::vector<int32_t> s_k, t;

        auto t0 = std::chrono::steady_clock::now();
        KeyGen(n, q, A, s_k, t);
        auto t1 = std::chrono::steady_clock::now();
        sum_keygen_us += std::chrono::duration_cast<
            std::chrono::microseconds>(t1 - t0).count();

        std::vector<std::vector<int32_t>> AT;
        transposeA(A, AT);

        // ---- Encaps --------------------------------------------------------
        std::vector<int32_t> c, K_enc;

        auto t2 = std::chrono::steady_clock::now();
        Encaps(n, q, t, c, A, AT, K_enc);
        auto t3 = std::chrono::steady_clock::now();
        sum_encaps_us += std::chrono::duration_cast<
            std::chrono::microseconds>(t3 - t2).count();

        // ---- Decaps --------------------------------------------------------
        std::vector<int32_t> K_dec;

        auto t4 = std::chrono::steady_clock::now();
        Decaps(n, q, t, s_k, c, K_dec, A, AT);
        auto t5 = std::chrono::steady_clock::now();
        sum_decaps_us += std::chrono::duration_cast<
            std::chrono::microseconds>(t5 - t4).count();

        // ---- Correctness check ---------------------------------------------
        bool same = (K_enc.size() == K_dec.size());
        if (same)
            for (size_t i = 0; i < K_enc.size() && same; ++i)
                same = (K_enc[i] == K_dec[i]);
        if (!same) ++mismatches;
    }

    const auto wall_end = std::chrono::steady_clock::now();
    const long long wall_us = std::chrono::duration_cast<
        std::chrono::microseconds>(wall_end - wall_start).count();

    const double total_s       = wall_us / 1e6;
    const double avg_keygen_us = static_cast<double>(sum_keygen_us) / N;
    const double avg_encaps_us = static_cast<double>(sum_encaps_us) / N;
    const double avg_decaps_us = static_cast<double>(sum_decaps_us) / N;

    // ---- CSV output (format unchanged from legacy) -------------------------
    std::cout << "Results (CSV format):\n"
              << "n;threads;avg_keygen_us;avg_encaps_us;avg_decaps_us;"
                 "total_s;mismatches\n"
              << n                     << ";"
              << omp_get_max_threads() << ";"
              << avg_keygen_us         << ";"
              << avg_encaps_us         << ";"
              << avg_decaps_us         << ";"
              << total_s               << ";"
              << mismatches            << "\n";

    // ---- Mismatch warning and suggestions ----------------------------------
    if (mismatches > 0 && cfg.warn_mismatch) {
        std::cerr << "\nWARNING: " << mismatches << " / " << N
                  << " key mismatch(es) detected.\n";

        if (cfg.suggest_fix) {
            std::cerr
                << "\nSUGGESTION: The most common cause for n=" << n
                << ", q=" << q << " is noise that is too large.\n"
                << "Add or update the following in your config file:\n"
                << "  noise_kind  = binomial_eta1   # sets e, r, e1, e2 all to eta=1\n"
                << "  s_eta       = 3               # secret key eta (usually keep at 3)\n"
                << "\nFor a custom investigation, set individual noise kinds:\n"
                << "  e_kind   = binomial_eta1\n"
                << "  r_kind   = binomial_eta1\n"
                << "  e1_kind  = binomial_eta1\n"
                << "  e2_kind  = binomial_eta1\n"
                << "\nSee validate_config output above for the exact noise margin.\n";
        }
    }

    return (mismatches == 0) ? 0 : 1;
}
