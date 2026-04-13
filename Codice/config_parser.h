#pragma once
// config_parser.h
// Declarations for the config file parser and parameter validator.
// ============================================================================

#include "config.h"
#include <string>

// ---------------------------------------------------------------------------
// parse_config
// ---------------------------------------------------------------------------
// Reads a plain-text key=value config file (see sample config.txt).
// Lines starting with '#' or empty lines are ignored.
// Keys are case-insensitive; unknown keys emit a warning to stderr.
//
// Returns a fully populated LweKemConfig.
// Throws std::runtime_error on unrecoverable parse errors.
LweKemConfig parse_config(const std::string &path);

// ---------------------------------------------------------------------------
// parse_cli_fallback
// ---------------------------------------------------------------------------
// Provides backward-compatible CLI parsing:  ./lwe_kem N n
// Returns a default LweKemConfig with N and n filled from argv.
LweKemConfig parse_cli_fallback(int argc, char **argv);

// ---------------------------------------------------------------------------
// validate_config
// ---------------------------------------------------------------------------
// Checks parameter constraints and prints warnings.
// Returns false (and prints a detailed explanation) if the noise parameters
// are incompatible with (n, q) — i.e. decryption failures are expected.
// Does NOT abort; the caller decides whether to continue.
bool validate_config(const LweKemConfig &cfg, std::ostream &out);

// ---------------------------------------------------------------------------
// print_config
// ---------------------------------------------------------------------------
// Pretty-prints all parameters to `out` for reproducibility headers.
void print_config(const LweKemConfig &cfg, std::ostream &out);
