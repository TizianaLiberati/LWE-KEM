#include "openfhe.h"
#include "openfhecore.h"

#include "lwe-pke.h"
#include "binfhecontext.h"

#include "utils/hashutil.h"

#include "math/discretegaussiangenerator.h"
#include "math/distributiongenerator.h"
#include "math/hal/integer.h"
#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>

#include <random>
#include <sstream>
#include <iomanip>

#include <chrono>

#include <bitset>

#include <omp.h>

using namespace lbcrypto;

int mod(int a, int b)
{
    return (a % b + b) % b;
}

std::vector<std::vector<int32_t>> GenerateRandomMatrixInt32(size_t n, size_t m, int32_t maxValue)
{
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int32_t> dist(0, maxValue);

    std::vector<std::vector<int32_t>> matrix(n, std::vector<int32_t>(m));

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            matrix[i][j] = dist(rng);
        }
    }

    return matrix;
}

std::vector<uint32_t> GenerateRandomBitVectorUInt32(size_t n)
{
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<uint32_t> dist(0, 1);

    std::vector<uint32_t> vec;
    vec.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        vec.push_back(dist(rng));
    }

    return vec;
}

std::vector<int32_t> GenerateGaussianVector(size_t m, NativeInteger q, double stddev)
{
    DiscreteGaussianGeneratorImpl<NativeVector> dgg(stddev);
    std::vector<int32_t> vec(m);
    int32_t q_int = static_cast<int32_t>(q.ConvertToInt()); // q in forma intera

    for (size_t i = 0; i < m; ++i)
    {
        auto sample_big = dgg.GenerateInteger(0, stddev, m);
        int32_t sample = sample_big;
        vec[i] = sample;
    }

    return vec;
}

std::vector<std::vector<int32_t>> transpose(const std::vector<std::vector<int32_t>> &matrix)
{
    int32_t rows = matrix.size();
    int32_t cols = matrix[0].size();
    std::vector<std::vector<int32_t>> result(cols, std::vector<int32_t>(rows));

    for (uint32_t i = 0; i < rows; ++i)
        for (uint32_t j = 0; j < cols; ++j)
            result[j][i] = matrix[i][j];

    return result;
}

int32_t sample_eta_centered_binomial(uint8_t eta, std::mt19937 &gen)
{
    std::uniform_int_distribution<uint8_t> dis(0, 1); // bit 0 o 1
    uint8_t sum1 = 0, sum2 = 0;

    for (uint8_t i = 0; i < eta; ++i)
    {
        sum1 += dis(gen);
        sum2 += dis(gen);
    }

    return static_cast<int32_t>(sum1) - static_cast<int32_t>(sum2);
}

std::vector<int32_t> sample_vector_binomial(uint32_t n, uint8_t eta)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int32_t> result;
    result.reserve(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        result.push_back(sample_eta_centered_binomial(eta, gen));
    }

    return result;
}

int32_t getRandomInt(int min, int max)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(min, max);
    return distr(gen);
}

void KeyGen(uint32_t n, uint32_t m, uint32_t q, double stddev, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &s, std::vector<int32_t> &t, int32_t bound)
{
    A = GenerateRandomMatrixInt32(n, m, q - 1);
    s = sample_vector_binomial(n, bound);
    std::vector<int32_t> e = GenerateGaussianVector(m, q, stddev);
    std::vector<int32_t> prod(m, q);

    for (uint32_t i = 0; i < m; ++i)
    {
        prod[i] = 0;
        for (uint32_t j = 0; j < n; ++j)
        {

            prod[i] = mod(prod[i] + A[i][j] * s[j], q);
        }
    }
    std::vector<int32_t> t1(m, 0);

    for (uint32_t i = 0; i < m; ++i)
    {
        t1[i] = mod(prod[i] + e[i], q);
    }

    t = t1;
}

void Encrypt(uint32_t n, uint32_t m, uint32_t q, double stddev, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2)
{
    std::vector<int32_t> prod(m, q);
    std::vector<std::vector<int32_t>> A_transp = transpose(A);

    for (uint32_t i = 0; i < n; ++i)
    {
        prod[i] = 0;
        for (uint32_t j = 0; j < m; ++j)
        {
            prod[i] = mod(prod[i] + A_transp[i][j] * r[j], q);
        }
    }

    std::vector<int32_t> u1(n, q);
    for (uint32_t i = 0; i < n; ++i)
    {
        u1[i] = prod[i] + e1[i];
    }
    u = u1;

    int32_t v1 = 0;
    int32_t risultato = 0;
    for (size_t i = 0; i < t.size(); ++i)
    {
        risultato = mod(risultato + t[i] * r[i], q);
    }
    v1 = risultato + e2 + plaintext_i;
    v_i = v1;
}

void Decrypt(int32_t &v_i, std::vector<int32_t> u, std::vector<int32_t> s, uint32_t q, int32_t &decrypt_i, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2)
{
    int32_t risultato = 0;
    for (size_t i = 0; i < s.size(); ++i)
    {
        risultato = mod(risultato + (s[i] * u[i]), q);
    }
    int32_t mu = mod(v_i - risultato, q);

    uint32_t m;

    int32_t bound = q / 4;

    if (mu > (q - bound) || mu <= bound)
    {
        m = 0;
    }
    else
    {
        m = q / 2;
    }

    decrypt_i = m;
}

std::string Int32VectorToString(const std::vector<int32_t> &vec)
{
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i)
    {
        oss << vec[i];
        if (i != vec.size() - 1)
        {
            oss << ",";
        }
    }
    return oss.str();
}

std::string HashToBinaryString_256(const std::vector<int64_t> &digest)
{
    std::ostringstream binaryStream;
    for (auto val : digest)
    {
        uint8_t byte = static_cast<uint8_t>(val);
        binaryStream << std::bitset<8>(byte);
    }

    std::string binary = binaryStream.str();
    return binary.substr(0, 256);
}

std::vector<int32_t> BitStringToInt32Vector(const std::string &bitstring)
{
    std::vector<int32_t> result;
    if (mod(bitstring.size(), 8) != 0)
    {
        throw std::invalid_argument("La bitstring deve avere una lunghezza multipla di 8.");
    }

    for (size_t i = 0; i + 8 <= bitstring.size(); i += 8)
    {
        std::string byte_str = bitstring.substr(i, 8);
        uint32_t val = std::stoul(byte_str, nullptr, 2); // Conversione binaria
        result.push_back(static_cast<int32_t>(val));
    }

    return result;
}

std::vector<int32_t> FromMatrixToVector(const std::vector<std::vector<int32_t>> &matrix, std::vector<int32_t> &t)
{
    std::vector<int32_t> vec((matrix.size() * matrix[0].size()) + t.size(), 0);
    int count = 0;
    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix[i].size(); ++j)
        {
            vec[count] = matrix[i][j];
            count++;
        }
    }
    for (int i = 0; i < t.size(); ++i)
    {
        vec[matrix.size() * matrix[0].size() + i] = t[i];
    }
    return vec;
}

void Encaps(uint32_t n, uint32_t m, uint32_t q, double stddev, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, int32_t plaintext_i, std::vector<int32_t> &Hash_K, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2)
{
    std::vector<int64_t> digest_conc;
    std::vector<int64_t> digest_final;

    std::vector<int32_t> conc = FromMatrixToVector(A, t);

    // Hash di t
    std::string conc_string = Int32VectorToString(conc);
    HashUtil::Hash(conc_string, HashAlgorithm::SHA_256, digest_conc);

    std::string bitstring_conc = HashToBinaryString_256(digest_conc);

    std::vector<int32_t> Hash_conc = BitStringToInt32Vector(bitstring_conc);
    std::vector<int32_t> Hash_concm(Hash_conc.size() + 1, 0);

    for (size_t i = 0; i < Hash_conc.size(); ++i)
    {
        Hash_concm[i] = Hash_conc[i];
    }
    Hash_concm[Hash_conc.size()] = plaintext_i;

    std::string final_string = Int32VectorToString(Hash_concm);
    HashUtil::Hash(final_string, HashAlgorithm::SHA_256, digest_final);
    std::string bitstring_final = HashToBinaryString_256(digest_final);

    std::vector<int32_t> K_caps = BitStringToInt32Vector(bitstring_final);

    Encrypt(n, m, q, stddev, A, t, u, v_i, plaintext_i, r, e1, e2);

    std::vector<int32_t> c(u.size() + 1, 0);
    for (int i = 0; i < u.size(); ++i)
    {
        c[i] = u[i];
    }
    c[u.size()] = v_i;

    std::vector<int64_t> digest_c;
    std::string c_string = Int32VectorToString(c);
    HashUtil::Hash(c_string, HashAlgorithm::SHA_256, digest_c);
    std::string bitstring_c = HashToBinaryString_256(digest_c);
    std::vector<int32_t> Hash_c = BitStringToInt32Vector(bitstring_c);

    std::vector<int32_t> K(K_caps.size() + Hash_c.size(), 0);
    for (size_t i = 0; i < K_caps.size(); ++i)
    {
        K[i] = K_caps[i];
    }
    for (int i = 0; i < Hash_c.size(); ++i)
    {
        K[K_caps.size() + i] = Hash_c[i];
    }

    std::vector<int64_t> digest_K;
    std::string K_string = Int32VectorToString(K);
    HashUtil::Hash(K_string, HashAlgorithm::SHA_256, digest_K);
    std::string bitstring_K = HashToBinaryString_256(digest_K);
    Hash_K = BitStringToInt32Vector(bitstring_K);
}

void Decaps(int32_t &v_i, std::vector<int32_t> &u, std::vector<int32_t> s, uint32_t q, int32_t &decrypt_i, std::vector<int32_t> &t, uint32_t n, uint32_t m, double stddev, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2)
{
    uint32_t m_dec;
    Decrypt(v_i, u, s, q, decrypt_i, r, e1, e2);
    m_dec = decrypt_i;

    std::vector<int64_t> digest_conc;
    std::vector<int64_t> digest_final;

    std::vector<int32_t> conc = FromMatrixToVector(A, t);

    // Hash di t
    std::string conc_string = Int32VectorToString(conc);
    HashUtil::Hash(conc_string, HashAlgorithm::SHA_256, digest_conc);
    std::string bitstring_conc = HashToBinaryString_256(digest_conc);

    std::vector<int32_t> Hash_conc = BitStringToInt32Vector(bitstring_conc);
    std::vector<int32_t> Hash_concm(Hash_conc.size() + 1, 0);

    for (size_t i = 0; i < Hash_conc.size(); ++i)
    {
        Hash_concm[i] = Hash_conc[i];
    }
    Hash_concm[Hash_conc.size()] = m_dec;

    std::string final_string = Int32VectorToString(Hash_concm);
    HashUtil::Hash(final_string, HashAlgorithm::SHA_256, digest_final);
    std::string bitstring_final = HashToBinaryString_256(digest_final);

    std::vector<int32_t> K_caps = BitStringToInt32Vector(bitstring_final);

    std::vector<int32_t> u_new(u.size(), 0);
    int32_t v_i_new = 0;
    Encrypt(n, m, q, stddev, A, t, u_new, v_i_new, m_dec, r, e1, e2);

    std::vector<int32_t> c(u_new.size() + 1, 0);
    for (int i = 0; i < u_new.size(); ++i)
    {
        c[i] = u_new[i];
    }
    c[u_new.size()] = v_i_new;

    std::vector<int64_t> digest_c;
    std::string c_string = Int32VectorToString(c);
    HashUtil::Hash(c_string, HashAlgorithm::SHA_256, digest_c);
    std::string bitstring_c = HashToBinaryString_256(digest_c);
    std::vector<int32_t> Hash_c = BitStringToInt32Vector(bitstring_c);

    std::vector<int32_t> K(K_caps.size() + Hash_c.size(), 0);
    for (size_t i = 0; i < K_caps.size(); ++i)
    {
        K[i] = K_caps[i];
    }
    for (int i = 0; i < Hash_c.size(); ++i)
    {
        K[K_caps.size() + i] = Hash_c[i];
    }

    std::vector<int64_t> digest_K;
    std::string K_string = Int32VectorToString(K);
    HashUtil::Hash(K_string, HashAlgorithm::SHA_256, digest_K);
    std::string bitstring_K = HashToBinaryString_256(digest_K);
    auto Hash_K = BitStringToInt32Vector(bitstring_K);
    /* Utile per testare se la decaps funziona
    if (v_i_new == v_i && u_new == u)
    {
        std::cout << "Decapsulation success, ss: " << Hash_K << std::endl;
    }
    else
    {
        std::cout << "Decapsulation Failed" << std::endl;
    }
    */
}

int main()
{
    uint32_t n = 1024; // Da paper
    uint32_t m = 1024; // Da paper

    double stddev = 2.3; // Da paper

    int bound = 3;

    uint32_t q = 3329; // Da paper

    double GlobalTime;

    std::vector<uint32_t> plaintext = GenerateRandomBitVectorUInt32(n);
    for (uint32_t i = 0; i < plaintext.size(); ++i)
    {
        if (plaintext[i] == 0)
        {
            plaintext[i] = 0;
        }
        else
        {
            plaintext[i] = q / 2;
        };
    }

    uint32_t plaintext_i;

    std::vector<uint32_t> decrypt(plaintext.size(), 0);
    int32_t decrypt_i = 0;

    std::vector<int32_t> s(n, 0);

    std::vector<int32_t> t(m, 0);

    std::vector<int32_t> v(plaintext.size(), q);
    int32_t v_i;

    std::vector<int32_t> Hash_K;
    std::vector<std::vector<int32_t>> A;

    auto start = std::chrono::high_resolution_clock::now();
    KeyGen(n, m, q, stddev, A, s, t, bound);

    std::vector<int32_t> prova_v(plaintext.size(), 0);
#pragma omp parallel for num_threads(60)
    for (uint32_t i = 0; i < plaintext.size(); ++i)
    {

        std::vector<int32_t> r = GenerateGaussianVector(n, q, stddev);
        std::vector<int32_t> e1 = GenerateGaussianVector(n, q, stddev);

        int32_t e2 = getRandomInt((-1) * bound, bound);

        plaintext_i = plaintext[i];
        std::vector<int32_t> u(n, q);
        Encaps(n, m, q, stddev, A, t, u, v_i, plaintext_i, Hash_K, r, e1, e2);
        Decaps(v_i, u, s, q, decrypt_i, t, n, m, stddev, A, r, e1, e2);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto GlobalTimeFor = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << n << ";" << stddev << ";" << bound << ";" << q << ";" << GlobalTimeFor << std::endl;

    return 0;
}
