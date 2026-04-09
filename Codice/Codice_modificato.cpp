#include <array>
#include <random>
#include <limits>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>

#include <cstdlib> //NEW per makefile automatizzato

#include <omp.h> // openMP

#include "pke.h" // funzioni del pke
#include "utils.h" // funzioni ausiliarie
#include "hash_openssl.h" // stiamo usando questa
#include "noise.h" // funzioni di generazione dei noise
#include "kem.h" // funzioni del kem



int main(int argc, char** argv)
{
    uint32_t q = 3329;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <N> <n>\n";
        return 1;
    }
    int N = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    size_t sizeN = N * sizeof(double);
    size_t sizen = n * sizeof(double);
    std::cout << "N = " << N << ", sizeN = " << sizeN << std::endl;
    std::cout << "n = " << n << ", sizen = " << sizen << std::endl;

    auto startTot = std::chrono::steady_clock::now();
    
    long long sum_keygen_us = 0;
    long long sum_encaps_us = 0;
    long long sum_decaps_us = 0;

    int mismatches = 0;

    for (int k = 0; k < N; ++k){
        
        std::vector<std::vector<int32_t>> A;
        std::vector<int32_t> s_k, t;

        auto startK = std::chrono::steady_clock::now();
        KeyGen(n, q, A, s_k, t);
        auto endK = std::chrono::steady_clock::now();

        auto elapsedK = std::chrono::duration_cast<std::chrono::microseconds>(endK - startK);
        sum_keygen_us += elapsedK.count();        

        std::vector<std::vector<int32_t>> AT;
        transposeA(A, AT);

        std::vector<int32_t> c;
        std::vector<int32_t> K_enc; // chiave condivisa generata da utente che cifra
    
        auto startE = std::chrono::steady_clock::now();
        Encaps(n, q, t, c, A, AT, K_enc); // K_enc ha 32 byte = 256 bit
        auto endE = std::chrono::steady_clock::now();
        auto elapsedE = std::chrono::duration_cast<std::chrono::microseconds>(endE - startE);
        sum_encaps_us += elapsedE.count();

        std::vector<int32_t> K_dec; // chiave condivisa generata da utente che decifra

        auto startD = std::chrono::steady_clock::now();
        Decaps(n, q, t, s_k, c, K_dec, A, AT);
        auto endD = std::chrono::steady_clock::now();
        auto elapsedD = std::chrono::duration_cast<std::chrono::microseconds>(endD - startD);
        sum_decaps_us += elapsedD.count();
        
        bool sameK = (K_enc.size() == K_dec.size());
        if (sameK)
        {
            for (size_t i = 0; i < K_enc.size(); ++i)
            {
                if (K_enc[i] != K_dec[i])
                {
                    sameK = false;
                    break;
                }
            }
            if (!sameK) mismatches++;
        }
       
    }

    auto endTot = std::chrono::steady_clock::now();

    // microsecondi = 10^(-6) secondi
    auto elapsedTot_us = std::chrono::duration_cast<std::chrono::microseconds>(endTot - startTot).count();
    double total_s = elapsedTot_us / 1e6;

    double avg_keygen_us = (double)sum_keygen_us / N;
    double avg_encaps_us = (double)sum_encaps_us / N;
    double avg_decaps_us = (double)sum_decaps_us / N;
    
    int threads = omp_get_max_threads();
    std::cout << n << ";" << threads << ";"
              << avg_keygen_us << ";"
              << avg_encaps_us << ";"
              << avg_decaps_us << ";"
              << total_s 
              << "\n";
    
    return 0;
}