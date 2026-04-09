#pragma once

#include <cstdint>
#include <limits>

inline Botan::RandomNumberGenerator& global_rng()
{
    static Botan::AutoSeeded_RNG rng; //variabile statica
    return rng;
}

// Adattatore per usare global_rng() con <random>, rende compatibile
/*  
    prima con le funzioni usate per le distribuzioni prendevano in input un 
    UniformRandomBitGenerator (URBG), per non cambiare quelle funzioni definiamo 
    la seguente struttura che rispetta URBG
*/
struct Botan_URBG
{
    using result_type = std::uint32_t;

    Botan_URBG()
        : rng_(global_rng())
    {}

    static constexpr result_type min()
    {
        return std::numeric_limits<result_type>::min();
    }

    static constexpr result_type max()
    {
        return std::numeric_limits<result_type>::max();
    }

    result_type operator()()
    {
        result_type x;
        rng_.randomize(reinterpret_cast<std::uint8_t*>(&x), sizeof(x));
        return x;
    }

private:
    Botan::RandomNumberGenerator& rng_;
};