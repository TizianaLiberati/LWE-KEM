#include <vector>
#include <cstdint>

// #include <botan/system_rng.h> 

#include "kem.h"
#include "utils.h"
#include "hash.h"
#include "noise.h"
#include "pke.h"

/////////////////////////////////////   Encaps  /////////////////////////////////////

void Encaps(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &c, const std::vector<std::vector<int32_t>> &A, const std::vector<std::vector<int32_t>> &AT, std::vector<int32_t> &Hash_K)
{
    const size_t msg_bits = 256;

    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh = SHA3_256(At); // potrei far in modo che SHA3_256 restituisca in uint8_t (magari se passiamo in SHA3-526)

    // Genero il plaintext lungo 256 random
    std::vector<int32_t> m(msg_bits, 0);
    // std::mt19937 rng(std::random_device{}());



    // TODO: sostituisci mt19937 con generatore Botan o con classe <random> 
    // Botan::System_RNG rng;
    
    // std::uniform_int_distribution<int> bit(0, 1);
    // for (uint32_t i = 0; i < msg_bits; ++i)
    //      m[i] = bit(rng);

    for (uint32_t i = 0; i < msg_bits; ++i)
        m[i] = random_bit(); //nuova funzione definita in utils.cpp che campiona randomicamente un bit

    // Creo il seed che userò nella XOF
    std::vector<int32_t> pkh_m = concat(pkh, m);
    std::vector<int32_t> seed = SHA3_256(pkh_m);


    // implementa XOF e prendi i seed per i noise
    std::vector<uint8_t> coins = xof_coins(seed, n, msg_bits);


    // Li genero dopo 

    // std::vector<int32_t> r = GenerateGaussianVector(n);
    // std::vector<int32_t> e1 = GenerateGaussianVector(n);

    // int32_t e2 = getRandomInt((-1) * 3, 3); 

    std::vector<int32_t> K_cap = SHA3_256(concat(pkh, m));

    c.assign((size_t)msg_bits * (n + 1), 0);

    std::vector<int32_t> u(n, 0);
    int32_t v_i = 0;

    size_t pos = 0; // indice dentro coins

    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        // genera r, e1, e2 per ogni bit del plaintext (nel for)
        NoiseTriple noise_i = GenerateNoisesForOneBit(coins, pos, n);

        std::vector<int32_t>& r  = noise_i.r;
        std::vector<int32_t>& e1 = noise_i.e1;
        int32_t e2 = noise_i.e2;

        int32_t plaintext_j = (m[j] == 1) ? (int32_t)(q / 2) : 0;

        Encrypt(n, q, t, u, v_i, plaintext_j, r, e1, e2, AT);
        size_t off = (size_t)(j) * (n + 1);

        for (uint32_t i = 0; i < n; ++i)
        {
            c[off + i] = u[i];
        }

        c[off + n] = v_i;
    }

    std::vector<int32_t> h_c = SHA3_256(c);
    Hash_K = SHA3_256(concat(K_cap, h_c));
}

/////////////////////////////////////   Decaps  /////////////////////////////////////

void Decaps(uint32_t n, uint32_t q, const std::vector<int32_t> &t, const std::vector<int32_t> &s_k, const std::vector<int32_t> &c, std::vector<int32_t> &Hash_K, const std::vector<std::vector<int32_t>> &A, const std::vector<std::vector<int32_t>> &AT)
{
    const size_t msg_bits = 256;

    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh = SHA3_256(At);

    std::vector<int32_t> mprime(msg_bits, 0);
    std::vector<int32_t> u_j(n, 0);
    int32_t v_j = 0, dec_m = 0;

    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        size_t off = (size_t)j * (n + 1);

        for (uint32_t i = 0; i < n; ++i)
            u_j[i] = c[off + i];

        v_j = c[off + n];
        Decrypt(v_j, u_j, s_k, q, dec_m);
        mprime[j] = (dec_m == (int32_t)(q / 2)) ? 1 : 0;
    }

    std::vector<int32_t> K_cap = SHA3_256(concat(pkh, mprime));

    std::vector<int32_t> cchk;
    cchk.assign((size_t)msg_bits * (n + 1), 0);
    std::vector<int32_t> u_tmp(n, 0);
    int32_t v_tmp = 0;

    // Creo il seed che userò nella XOF
    std::vector<int32_t> pkh_mprime = concat(pkh, mprime);
    std::vector<int32_t> seed = SHA3_256(pkh_mprime);

    // implementa XOF e prendi i seed per i noise
    std::vector<uint8_t> coins = xof_coins(seed, n, msg_bits);
    size_t pos = 0; // indice dentro coins

    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        int32_t m_j_map = mprime[j] ? (int32_t)(q / 2) : 0;

        // TODO: ricontrolla perchè li generi dopo

        // std::vector<int32_t> r = GenerateGaussianVector(n);
        // std::vector<int32_t> e1 = GenerateGaussianVector(n);
        // int32_t e2 = getRandomInt(-3, 3);

        // genera r, e1, e2 per ogni bit del plaintext (nel for)
        NoiseTriple noise_i = GenerateNoisesForOneBit(coins, pos, n);

        std::vector<int32_t>& r  = noise_i.r;
        std::vector<int32_t>& e1 = noise_i.e1;
        int32_t e2 = noise_i.e2;
        Encrypt(n, q, const_cast<std::vector<int32_t> &>(t), u_tmp, v_tmp, (uint32_t)m_j_map, r, e1, e2, 
                const_cast<std::vector<std::vector<int32_t>> &>(AT));

        size_t off = (size_t)j * (n + 1);
        for (uint32_t i = 0; i < n; ++i)
            cchk[off + i] = u_tmp[i];
        cchk[off + n] = v_tmp;
    }

    bool equal = (cchk.size() == c.size());
    if (equal)
    {
        for (size_t k = 0; k < c.size(); ++k)
            if (cchk[k] != c[k])
            {
                equal = false;
                break;
            }
    }
    
    std::vector<int32_t> h_c = SHA3_256(c);
    if (equal)
    {
        Hash_K = SHA3_256(concat(K_cap, h_c));
        // std::cout << "ha funzionato\n";
    }
    else
    {
        // implicit rejection
        std::vector<int32_t> z(s_k.begin() + n, s_k.end());
        Hash_K = SHA3_256(concat(z, h_c));
        // std::cout << "non ha funzionato\n";
    }
}