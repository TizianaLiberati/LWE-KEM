#include <vector>
#include <cstdint>
#include <iostream>

#include <omp.h> // openMP

#include "pke.h"
#include "utils.h"

/////////////////////////////////////   KeyGen  /////////////////////////////////////

void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &s_k, std::vector<int32_t> &t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, q);

    std::vector<int32_t> z(256); //256 bits

    #pragma omp parallel for //OKKKS
    for (int i = 0; i < 256; ++i)
        z[i] = getRandomInt(0, 1);

    // secret key s_k data dalla concat di s e z
    s_k = concat(s, z);

    // auto start_prod = std::chrono::steady_clock::now();
        
    // #pragma omp parallel for collapse(2) //funziona male con questo
    #pragma omp parallel for //OKKKS
    for (uint32_t i = 0; i < n; ++i)
    {
        prod[i] = 0;
        for (uint32_t j = 0; j < n; ++j)
        {

            prod[i] = mod(prod[i] + A[i][j] * s[j], q);
        }
    }
    // auto end_prod = std::chrono::steady_clock::now();
    // auto elapsed_prod = std::chrono::duration_cast<std::chrono::microseconds>(end_prod - start_prod);
    // std::cout << "Prod time: " << elapsed_prod.count() << " mus\n";

    std::vector<int32_t> t1(n, 0);
    #pragma omp parallel for //OKKKS
    for (uint32_t i = 0; i < n; ++i)
    {
        t1[i] = mod(prod[i] + e[i], q);
    }
    
    t = t1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Encrypt /////////////////////////////////////
void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2, const std::vector<std::vector<int32_t>> &AT)
{
    std::vector<int32_t> prod(n, q);

    // #pragma omp parallel for collapse(2)
    #pragma omp parallel for //OKKKS
    for (uint32_t i = 0; i < n; ++i)
    {
        prod[i] = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            prod[i] = mod(prod[i] + AT[i][j] * r[j], q);
        }
    }

    std::vector<int32_t> u1(n, q);
    #pragma omp parallel for //OKKKS
    for (uint32_t i = 0; i < n; ++i)
    {
        u1[i] = mod(prod[i] + e1[i], q);
    }
    u = u1;

    int32_t v1 = 0;
    int32_t risultato = 0;
    
    for (size_t i = 0; i < t.size(); ++i)
    {
        risultato = mod(risultato + t[i] * r[i], q);
    }
    v1 = mod(risultato + e2 + plaintext_i, q);
    v_i = v1;
}
/////////////////////////////////////   Decrypt /////////////////////////////////////

void Decrypt(int32_t v_i, const std::vector<int32_t> &u, const std::vector<int32_t> &s_k, uint32_t q, int32_t &decrypt_i)
{
    long long dot = 0;

    // Considero solo la prima parte di s_k, quella relativa ad s
    const size_t n = u.size(); // la dimensione di s (senza z) è n = dimensione di u
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n); // definisco s uguale al primo elemento di s_k fino all'elemento n-esimo di s_k

    #pragma omp parallel for reduction(+:dot)  //OKKKS
    for (size_t i = 0; i < s.size(); ++i)
        dot += (long long)s[i] * (long long)u[i];

    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0)
        r += (long long)q;
    int32_t mu = (int32_t)r;

    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)q / 2;
}