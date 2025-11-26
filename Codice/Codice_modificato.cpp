#include <array>
#include <random>
#include <limits>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>

#include "pke.h" // funzioni del pke
#include "utils.h" // funzioni ausiliarie
#include "hash.h" // funzioni di hash
#include "noise.h" // funzioni di generazione dei noise
#include "kem.h" // funzioni del kem


int main()
{
    auto startTot = std::chrono::steady_clock::now();

    uint32_t n = 512;
    uint32_t q = 3329;

    std::vector<std::vector<int32_t>> A;
    std::vector<int32_t> s_k, t;

    auto startK = std::chrono::steady_clock::now();
    KeyGen(n, q, A, s_k, t);
    auto endK = std::chrono::steady_clock::now();

    auto elapsedK = std::chrono::duration_cast<std::chrono::microseconds>(endK - startK);
    std::cout << "KeyGen time: " << elapsedK.count() << " mus\n";

    std::vector<std::vector<int32_t>> AT;
    transposeA(A, AT);

    std::vector<int32_t> c;
    std::vector<int32_t> K_enc; // chiave condivisa generata da utente che cifra
  
    auto startE = std::chrono::steady_clock::now();
    Encaps(n, q, t, c, A, AT, K_enc); // K_enc ha 32 byte = 256 bit
    auto endE = std::chrono::steady_clock::now();
    auto elapsedE = std::chrono::duration_cast<std::chrono::microseconds>(endE - startE);
    std::cout << "Encaps time: " << elapsedE.count() << " mus\n";

    std::vector<int32_t> K_dec; // chiave condivisa generata da utente che decifra

    auto startD = std::chrono::steady_clock::now();
    Decaps(n, q, t, s_k, c, K_dec, A, AT);
    auto endD = std::chrono::steady_clock::now();
    auto elapsedD = std::chrono::duration_cast<std::chrono::microseconds>(endD - startD);
    std::cout << "Decaps time: " << elapsedD.count() << " mus\n";

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

    std::cout << "K_enc = [ ";
    for (auto x : K_enc)
        std::cout << x << " ";
    std::cout << "]\n";

    std::cout << "K_dec = [ ";
    for (auto x : K_dec)
        std::cout << x << " ";
    std::cout << "]\n";

    auto endTot = std::chrono::steady_clock::now();

    auto elapsedTot = std::chrono::duration_cast<std::chrono::microseconds>(endTot - startTot);
    std::cout << "Tot time: " << elapsedTot.count() << " mus\n";
    return 0;
}