// vedi anche https://github.com/Duthomhas/CSPRNG repo github

/* TODO:
    - vedi altra libreria per hash function: https://github.com/KeccakTeam/KeccakCodePackage (Team Keccak)
    https://github.com/KeccakTeam/KeccakCodePackage (openSSL)
    - arcora i seed
*/

#include <array>
#include <random>
#include <limits>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>

#include <omp.h> // openMP

#include "pke.h" // funzioni del pke
#include "utils.h" // funzioni ausiliarie
#include "hash.h" // funzioni di hash
#include "hash_openssl.h"
// #include "hash_keccak.h"
#include "noise.h" // funzioni di generazione dei noise
#include "kem.h" // funzioni del kem



int main()
{
    int N = 750;
    // auto startTot = std::chrono::steady_clock::now();

    uint32_t n = 512;
    uint32_t q = 3329;

    /*
    using clock_t = std::chrono::steady_clock; // è per il benchmark

    const size_t msg_bits = 256;

    // 1) input fisso come "stringa"
    std::string s = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";

    // 2) converti la stringa in vector<int32_t> (1 int32 per carattere)
    std::vector<int32_t> in;
    in.reserve(s.size());
    for (unsigned char ch : s)
        in.push_back((int32_t)ch);

    // warm-up (non misurato)
    for (int i = 0; i < 10; ++i) {
        auto seed = SHA3_256(in);
        auto coins = xof_coins(seed, n, msg_bits);
        (void)coins;
    }

    const int R = 1000; // numero prove
    long long sum_hash_us = 0;
    long long sum_xof_us  = 0;

    for (int i = 0; i < R; ++i) {
        auto t0 = clock_t::now();
        auto seed = SHA3_256(in);
        auto t1 = clock_t::now();
        auto coins = xof_coins(seed, n, msg_bits);
        auto t2 = clock_t::now();

        sum_hash_us += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        sum_xof_us  += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        if (coins.size() == 123456789) std::cout << "x"; // anti-ottimizzazione
    }

    std::cout << "[avg] SHA3_256 = " << (sum_hash_us / (double)R) << " ns\n";
    std::cout << "[avg] XOF(shake) = " << (sum_xof_us  / (double)R) << " us\n";

    // ---------------- BENCH OPENSSL ----------------

    // warm-up (non misurato)
    for (int i = 0; i < 10; ++i) {
        auto seed = SHA3_256_openssl(in);
        auto coins = xof_coins_openssl(seed, n, msg_bits);
        (void)coins;
    }

    long long sum_hash_ns_ossl = 0;
    long long sum_xof_us_ossl  = 0;

    for (int i = 0; i < R; ++i) {
        auto t0 = clock_t::now();
        auto seed = SHA3_256_openssl(in);
        auto t1 = clock_t::now();
        auto coins = xof_coins_openssl(seed, n, msg_bits);
        auto t2 = clock_t::now();

        sum_hash_ns_ossl += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        sum_xof_us_ossl  += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        if (coins.size() == 123456789) std::cout << "x";
    }

    std::cout << "[OpenSSL avg] SHA3_256 = " << (sum_hash_ns_ossl / (double)R) << " ns\n";
    std::cout << "[OpenSSL avg] XOF(shake) = " << (sum_xof_us_ossl  / (double)R) << " us\n";
*/

    auto startTot = std::chrono::steady_clock::now();
    #pragma omp parallel for
    

    for (int k = 0; k < N; ++k){
        
        std::vector<std::vector<int32_t>> A;
        std::vector<int32_t> s_k, t;

        //auto startK = std::chrono::steady_clock::now();
        KeyGen(n, q, A, s_k, t);
        // auto endK = std::chrono::steady_clock::now();

        // auto elapsedK = std::chrono::duration_cast<std::chrono::microseconds>(endK - startK);
        // std::cout << "KeyGen time: " << elapsedK.count() << " mus\n";

        std::vector<std::vector<int32_t>> AT;
        transposeA(A, AT);

        std::vector<int32_t> c;
        std::vector<int32_t> K_enc; // chiave condivisa generata da utente che cifra
    
        //auto startE = std::chrono::steady_clock::now();
        Encaps(n, q, t, c, A, AT, K_enc); // K_enc ha 32 byte = 256 bit
        // auto endE = std::chrono::steady_clock::now();
        // auto elapsedE = std::chrono::duration_cast<std::chrono::microseconds>(endE - startE);
        // std::cout << "Encaps time: " << elapsedE.count() << " mus\n";

        std::vector<int32_t> K_dec; // chiave condivisa generata da utente che decifra

        //auto startD = std::chrono::steady_clock::now();
        Decaps(n, q, t, s_k, c, K_dec, A, AT);
        // auto endD = std::chrono::steady_clock::now();
        // auto elapsedD = std::chrono::duration_cast<std::chrono::microseconds>(endD - startD);
        // std::cout << "Decaps time: " << elapsedD.count() << " mus\n";

        bool sameK = (K_enc.size() == K_dec.size());
        if (sameK)
        {
            for (size_t i = 0; i < K_enc.size(); ++i)
            {
                if (K_enc[i] != K_dec[i])
                {
                    std::cout << "Le chiavi sono diverse all'indice " << i << "\n";
                    sameK = false;
                    break;
                }
            }
        }

        // Stampe key encaps e decaps

        // std::cout << "K_enc = [ ";
        // for (auto x : K_enc)
        //     std::cout << x << " ";
        // std::cout << "]\n";

        // std::cout << "K_dec = [ ";
        // for (auto x : K_dec)
        //     std::cout << x << " ";
        // std::cout << "]\n";

         
    }

    auto endTot = std::chrono::steady_clock::now();

    // microsecondi = 10^(-6) secondi
    auto elapsedTot = std::chrono::duration_cast<std::chrono::seconds>(endTot - startTot);
    std::cout << "Tot time: " << elapsedTot.count() << " s\n";

    // auto endTot = std::chrono::steady_clock::now();

    // // microsecondi = 10^(-6) secondi
    // auto elapsedTot = std::chrono::duration_cast<std::chrono::seconds>(endTot - startTot);
    // std::cout << "Tot time: " << elapsedTot.count() << " s\n";

    
    return 0;
}