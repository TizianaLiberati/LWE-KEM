// config_parser.cpp
// Implementation of LWE-KEM config file parser, validator, and printer.
//
// File format (config.txt):
//   # comment
//   key = value
//
// Keys are case-insensitive.  Whitespace around '=' is ignored.
// Unknown keys produce a warning; they do not abort parsing.
// ============================================================================

#include "config_parser.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// ============================================================================
// Internal helpers
// ============================================================================

namespace {

// Trim leading and trailing whitespace from a string.
std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Convert string to lowercase.
std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

// Parse a NoiseKind from a string token.
// Accepted values (case-insensitive):
//   "binomial_eta1" | "eta1" | "1"
//   "binomial_eta3" | "eta3" | "3"
//   "gaussian"      | "gauss"
NoiseKind parse_noise_kind(const std::string &raw, const std::string &key) {
    std::string v = to_lower(trim(raw));
    if (v == "binomial_eta1" || v == "eta1" || v == "1")
        return NoiseKind::BINOMIAL_ETA1;
    if (v == "binomial_eta3" || v == "eta3" || v == "3")
        return NoiseKind::BINOMIAL_ETA3;
    if (v == "gaussian" || v == "gauss")
        return NoiseKind::GAUSSIAN;
    if (v == "uniform_centered_3" || v == "uniform_pm3" || v == "u3")
        return NoiseKind::UNIFORM_CENTERED_3;
    throw std::runtime_error("Unknown noise kind for key '" + key + "': " + raw
        + "\n  Valid values: binomial_eta1, binomial_eta3, gaussian");
}

// Compute sigma for a given NoiseKind/eta.
double noise_sigma(NoiseKind k, uint32_t eta, double gauss_sigma) {
    switch (k) {
        case NoiseKind::BINOMIAL_ETA1: return std::sqrt(0.5);
        case NoiseKind::BINOMIAL_ETA3: return std::sqrt(static_cast<double>(eta) / 2.0);
        case NoiseKind::GAUSSIAN:      return gauss_sigma;
        case NoiseKind::UNIFORM_CENTERED_3: return 2.0;
    }
    return 0.0;
}

} // anonymous namespace

// ============================================================================
// parse_config
// ============================================================================

LweKemConfig parse_config(const std::string &path) {
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open config file: " + path);

    LweKemConfig cfg;   // start from safe defaults
    int line_no = 0;

    std::string line;
    while (std::getline(f, line)) {
        ++line_no;
        std::string trimmed = trim(line);

        // Skip empty lines and comments.
        if (trimmed.empty() || trimmed[0] == '#') continue;

        // Split on first '='.
        size_t eq = trimmed.find('=');
        if (eq == std::string::npos) {
            std::cerr << "[config] Line " << line_no
                      << ": ignored (no '='): " << trimmed << "\n";
            continue;
        }

        std::string key = to_lower(trim(trimmed.substr(0, eq)));
        std::string val = trim(trimmed.substr(eq + 1));

        // Remove inline comments from value.
        size_t comment = val.find('#');
        if (comment != std::string::npos)
            val = trim(val.substr(0, comment));

        try {
            // ---- Problem parameters ----------------------------------------
            if      (key == "q")           cfg.q = static_cast<uint32_t>(std::stoul(val));
            else if (key == "n")           cfg.n = static_cast<uint32_t>(std::stoul(val));
            else if (key == "iterations" || key == "n_iterations" || key == "num_iterations" || key == "big_n")
                                           cfg.N = std::stoi(val);

            // ---- Noise parameters ------------------------------------------
            else if (key == "s_eta")       cfg.s_eta       = static_cast<uint32_t>(std::stoul(val));
            else if (key == "e_kind"  || key == "e_noise")  cfg.e_kind  = parse_noise_kind(val, key);
            else if (key == "r_kind"  || key == "r_noise")  cfg.r_kind  = parse_noise_kind(val, key);
            else if (key == "e1_kind" || key == "e1_noise") cfg.e1_kind = parse_noise_kind(val, key);
            else if (key == "e2_kind" || key == "e2_noise") cfg.e2_kind = parse_noise_kind(val, key);
            else if (key == "gauss_sigma" || key == "sigma") cfg.gauss_sigma = std::stod(val);

            // ---- Shortcut: set all encaps noise kinds at once ---------------
            else if (key == "noise_kind" || key == "all_noise") {
                NoiseKind k = parse_noise_kind(val, key);
                cfg.e_kind = cfg.r_kind = cfg.e1_kind = cfg.e2_kind = k;
            }

            // ---- Reproducibility -------------------------------------------
            else if (key == "seed")        cfg.seed = std::stoull(val);

            // ---- OpenMP ----------------------------------------------------
            else if (key == "num_threads" || key == "threads")
                                           cfg.num_threads = std::stoi(val);

            // ---- Output control --------------------------------------------
            else if (key == "verbose")     cfg.verbose       = (val == "1" || to_lower(val) == "true");
            else if (key == "warn_mismatch")  cfg.warn_mismatch = (val == "1" || to_lower(val) == "true");
            else if (key == "suggest_fix") cfg.suggest_fix   = (val == "1" || to_lower(val) == "true");

            // ---- Unknown key -----------------------------------------------
            else {
                std::cerr << "[config] Line " << line_no
                          << ": unknown key '" << key << "' — ignored.\n";
            }

        } catch (const std::exception &ex) {
            std::ostringstream oss;
            oss << "Parse error at line " << line_no
                << " (key='" << key << "', val='" << val << "'): " << ex.what();
            throw std::runtime_error(oss.str());
        }
    }

    return cfg;
}

// ============================================================================
// parse_cli_fallback
// ============================================================================

LweKemConfig parse_cli_fallback(int argc, char **argv) {
    LweKemConfig cfg;   // safe defaults

    if (argc >= 3) {
        cfg.N = std::atoi(argv[1]);
        cfg.n = static_cast<uint32_t>(std::atoi(argv[2]));
    }

    if (cfg.N <= 0) {
        throw std::runtime_error("N (iterations) must be a positive integer.");
    }
    if (cfg.n == 0) {
        throw std::runtime_error("n (lattice dimension) must be a positive integer.");
    }

    return cfg;
}

// ============================================================================
// validate_config
// ============================================================================

bool validate_config(const LweKemConfig &cfg, std::ostream &out) {
    bool ok = true;

    // Basic range checks.
    if (cfg.N <= 0) {
        out << "[WARN] iterations (N) must be > 0, got " << cfg.N << "\n";
        ok = false;
    }
    if (cfg.n == 0) {
        out << "[WARN] lattice dimension (n) must be > 0\n";
        ok = false;
    }
    if (cfg.q < 2) {
        out << "[WARN] modulus q must be >= 2\n";
        ok = false;
    }
    if (cfg.s_eta == 0 || cfg.s_eta > 5) {
        out << "[WARN] s_eta=" << cfg.s_eta << " is unusual; expected 1–3.\n";
    }

    // ---- Noise correctness check ------------------------------------------
    // Decryption noise term:  η = e·r + e2 − s·e1
    // For correct decryption: 6σ_η < q/4   (6-sigma safety margin)
    //
    // σ_η² = n·(σ_e² · σ_r²)  +  n·(σ_s² · σ_e1²)  +  σ_e2²
    // σ_s = sqrt(s_eta / 2)

    const double sigma_s  = std::sqrt(static_cast<double>(cfg.s_eta) / 2.0);
    const double sigma_e  = noise_sigma(cfg.e_kind,  1, cfg.gauss_sigma);
    const double sigma_r  = noise_sigma(cfg.r_kind,  1, cfg.gauss_sigma);
    const double sigma_e1 = noise_sigma(cfg.e1_kind, 1, cfg.gauss_sigma);
    const double sigma_e2 = noise_sigma(cfg.e2_kind, 1, cfg.gauss_sigma);

    const double var_eta = static_cast<double>(cfg.n) *
                           (sigma_e * sigma_e * sigma_r * sigma_r
                          + sigma_s * sigma_s * sigma_e1 * sigma_e1)
                         + sigma_e2 * sigma_e2;

    const double std_eta  = std::sqrt(var_eta);
    const double six_sigma = 6.0 * std_eta;
    const double q4        = static_cast<double>(cfg.q) / 4.0;

    out << std::fixed << std::setprecision(2);
    out << "[noise] σ_η = " << std_eta
        << "  6σ = " << six_sigma
        << "  q/4 = " << q4;

    if (six_sigma < q4) {
        out << "  ✓ SAFE (margin = " << (q4 - six_sigma) << ")\n";
    } else {
        out << "  ✗ UNSAFE — decryption failures EXPECTED\n";
        ok = false;

        out << "\n[MISMATCH LIKELY] With n=" << cfg.n
            << ", q=" << cfg.q << ":\n"
            << "  6σ_noise (" << six_sigma << ") >= q/4 (" << q4 << ")\n"
            << "  Per-bit failure prob ≈ "
            << 2.0 * std::erfc(q4 / std_eta / std::sqrt(2.0))
            << "\n";

        // Suggest safe alternatives.
        if (cfg.suggest_fix) {
            out << "\n[SUGGESTION] To fix mismatches for n=" << cfg.n
                << ", q=" << cfg.q << ":\n";

            // Option 1: reduce noise to eta=1 for all encaps terms.
            double sig1 = std::sqrt(0.5);   // eta=1
            double var1 = static_cast<double>(cfg.n) *
                          (sig1*sig1*sig1*sig1 + sigma_s*sigma_s*sig1*sig1)
                        + sig1*sig1;
            double six1 = 6.0 * std::sqrt(var1);
            out << "  Option A: set all noise to binomial_eta1 in config.txt:\n"
                << "            noise_kind = binomial_eta1\n"
                << "            → 6σ = " << six1 << (six1 < q4 ? " ✓" : " ✗") << "\n";

            // Option 2: increase q.
            uint32_t q_needed = static_cast<uint32_t>(
                std::ceil(4.0 * six_sigma)) + 1;
            out << "  Option B: increase q to at least " << q_needed
                << " in config.txt:\n"
                << "            q = " << q_needed << "\n";

            // Option 3: reduce n.
            // Find max n such that 6sigma < q/4 with current noise kinds.
            // var = n*(σ_e²σ_r² + σ_s²σ_e1²) + σ_e2²
            // n < ((q/4/6)² - σ_e2²) / (σ_e²σ_r² + σ_s²σ_e1²)
            double per_n = sigma_e*sigma_e*sigma_r*sigma_r
                         + sigma_s*sigma_s*sigma_e1*sigma_e1;
            if (per_n > 0) {
                double max_n = ((q4/6.0)*(q4/6.0) - sigma_e2*sigma_e2) / per_n;
                if (max_n > 0)
                    out << "  Option C: reduce n to at most "
                        << static_cast<uint32_t>(max_n)
                        << " in config.txt:\n"
                        << "            n = " << static_cast<uint32_t>(max_n) << "\n";
            }
        }
    }

    return ok;
}

// ============================================================================
// print_config
// ============================================================================

void print_config(const LweKemConfig &cfg, std::ostream &out) {
    out << "╔═══════════════════════════════════════════════╗\n"
        << "║          LWE-KEM Runtime Configuration        ║\n"
        << "╚═══════════════════════════════════════════════╝\n";

    auto row = [&](const char *label, auto val) {
        out << "  " << std::left << std::setw(22) << label << " = " << val << "\n";
    };

    row("modulus q",         cfg.q);
    row("dimension n",       cfg.n);
    row("iterations N",      cfg.N);
    out << "  ──────────────────────────────────────────────\n";
    row("secret eta (s)",    cfg.s_eta);
    row("keygen error (e)",  noise_kind_str(cfg.e_kind));
    row("encaps random (r)", noise_kind_str(cfg.r_kind));
    row("encaps error (e1)", noise_kind_str(cfg.e1_kind));
    row("encaps scalar(e2)", noise_kind_str(cfg.e2_kind));
    if (cfg.e_kind  == NoiseKind::GAUSSIAN ||
        cfg.r_kind  == NoiseKind::GAUSSIAN ||
        cfg.e1_kind == NoiseKind::GAUSSIAN ||
        cfg.e2_kind == NoiseKind::GAUSSIAN)
        row("gaussian sigma",    cfg.gauss_sigma);
    out << "  ──────────────────────────────────────────────\n";
    row("seed",              cfg.seed == 0 ? "random (entropy)" : std::to_string(cfg.seed));
    row("num_threads",       cfg.num_threads == 0 ? "OMP default" : std::to_string(cfg.num_threads));
    row("verbose",           cfg.verbose       ? "true" : "false");
    row("warn_mismatch",     cfg.warn_mismatch ? "true" : "false");
    row("suggest_fix",       cfg.suggest_fix   ? "true" : "false");
    out << "═══════════════════════════════════════════════\n";
}
