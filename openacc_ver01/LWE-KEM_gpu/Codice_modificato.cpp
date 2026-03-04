#include <vector>
#include <array>
#include <cstdint>
#include <random>
#include <limits>
#include <chrono>
#include <ctime>
#include <iostream>
#include "sha256.h"

using simple_sha256::sha256_bytes;
using simple_sha256::i32_to_le_bytes;


int mod(int a, int b)
{
    return (a % b + b) % b;
}

static std::vector<int32_t> SHA256(const std::vector<int32_t>& in) {
    std::vector<uint8_t> bytes; bytes.reserve(in.size()*4);
    for (int32_t x : in) {
        uint32_t ux = (uint32_t)x;
        bytes.push_back((uint8_t)(ux      & 0xFF));
        bytes.push_back((uint8_t)((ux>>8) & 0xFF));
        bytes.push_back((uint8_t)((ux>>16)& 0xFF));
        bytes.push_back((uint8_t)((ux>>24)& 0xFF));
    }
    auto digest = simple_sha256::sha256_bytes(bytes); // 32 byte
    std::vector<int32_t> out(32);
    for (size_t i=0;i<32;++i) out[i] = (int32_t)digest[i];
    return out;
}

std::vector<int32_t> concat(const std::vector<int32_t>& a, const std::vector<int32_t>& b) {
    std::vector<int32_t> out;
    out.reserve(a.size() + b.size());
    for (size_t i = 0; i < a.size(); ++i) out.push_back(a[i]);
    for (size_t i = 0; i < b.size(); ++i) out.push_back(b[i]);
    return out;
}

std::vector<int32_t> flatten_matrix(const std::vector<std::vector<int32_t>>& A) {
    std::vector<int32_t> out;
    if (A.empty()) return out;
    const size_t n = A.size();
    out.reserve(n * A[0].size());
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            out.push_back(A[i][j]);
    return out;
}

int32_t sample_discrete_gaussian(double sigma) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::normal_distribution<double> norm(0.0, sigma);

    double x = norm(gen);

    double bound = 6.0 * sigma;
    if (x >  bound) x =  bound;
    if (x < -bound) x = -bound;

    long long k = std::llround(x);

    if (k > std::numeric_limits<int32_t>::max()) k = std::numeric_limits<int32_t>::max();
    if (k < std::numeric_limits<int32_t>::min()) k = std::numeric_limits<int32_t>::min();
    return static_cast<int32_t>(k);
}

std::vector<int32_t> GenerateGaussianVector(size_t n) {
    std::vector<int32_t> vec;
    vec.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        vec.push_back(sample_discrete_gaussian(2.3));
    }
    return vec;
}

int32_t getRandomInt(int min, int max)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(min, max);
    return distr(gen);
}

std::vector<std::vector<int32_t>> GenerateRandomMatrixInt32(size_t n, int32_t maxValue)
{
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int32_t> dist(0, maxValue);

    std::vector<std::vector<int32_t>> matrix(n, std::vector<int32_t>(n));

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            matrix[i][j] = dist(rng);
        }
    }

    return matrix;
}

void transposeA(const std::vector<std::vector<int32_t>>& A,
               std::vector<std::vector<int32_t>>& AT) {
    const size_t n = A.size();
    AT.assign(n, std::vector<int32_t>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            AT[j][i] = A[i][j];
}

int32_t sample_eta_centered_binomial(uint8_t eta, std::mt19937 &gen) {
    std::uniform_int_distribution<uint8_t> dis(0, 1);
    uint8_t sum1 = 0, sum2 = 0;
    for (uint8_t i = 0; i < eta; ++i) {
        sum1 += dis(gen);
        sum2 += dis(gen);
    }
    return (int32_t)sum1 - (int32_t)sum2;
}

std::vector<int32_t> sample_vector_binomial(uint32_t n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<int32_t> result;
    result.reserve(n);
    for (uint32_t i = 0; i < n; ++i)
        result.push_back(sample_eta_centered_binomial(3, gen)); 
    return result;
}

void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &s, std::vector<int32_t> &t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, q);

    for (uint32_t i = 0; i < n; ++i)
    {
        prod[i] = 0;
        for (uint32_t j = 0; j < n; ++j)
        {

            prod[i] = mod(prod[i] + A[i][j] * s[j], q);
        }
    }

    std::vector<int32_t> t1(n, 0);

    for (uint32_t i = 0; i < n; ++i)
    {
        t1[i] = mod(prod[i] + e[i], q);
    }

    t = t1;
}

void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2, const std::vector<std::vector<int32_t>>& AT)
{
    std::vector<int32_t> prod(n, q);

    for (uint32_t i = 0; i < n; ++i)
    {
        prod[i] = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            prod[i] = mod(prod[i] + AT[i][j] * r[j], q);
        }
    }

    std::vector<int32_t> u1(n, q);
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

void Encaps(uint32_t n, uint32_t q, std::vector<int32_t> &t, const std::vector<int32_t> &plaintext, std::vector<int32_t> &c, std::vector<int32_t> &Hash_K, std::vector<int32_t> h_At, const std::vector<std::vector<int32_t>>& AT)
{
    std::vector<int32_t> r = GenerateGaussianVector(n);
    std::vector<int32_t> e1 = GenerateGaussianVector(n);
 
    int32_t e2 = getRandomInt((-1) * 3, 3);
           
    std::vector<int32_t> K_cap = SHA256(concat(h_At, plaintext));

    c.assign((size_t)n*(n+1), 0);

    std::vector<int32_t> u(n, 0);
    int32_t v_i = 0;

    for(uint32_t j = 0; j < n; ++j){

        int32_t plaintext_j = (plaintext[j] == 1) ? (int32_t)(q / 2) : 0;

        Encrypt(n, q, t, u, v_i, plaintext_j, r, e1, e2, AT);
        size_t off = (size_t)(j) * (n + 1);

        for (uint32_t i = 0; i < n; ++i) {
            c[off + i] = u[i];
        }
        
        c[off + n] = v_i;
    }
    
    std::vector<int32_t> h_c  = SHA256(c);  
    Hash_K = SHA256(concat(K_cap, h_c));
}

/*
void Decrypt(int32_t &v_i, std::vector<int32_t> u, std::vector<int32_t> s, uint32_t q, int32_t &decrypt_i)
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
*/

void Decrypt(int32_t v_i, const std::vector<int32_t>& u, const std::vector<int32_t>& s, uint32_t q, int32_t& decrypt_i) {
    long long dot = 0;
    for (size_t i = 0; i < s.size(); ++i)
        dot += (long long)s[i] * (long long)u[i];

    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0) r += (long long)q;
    int32_t mu = (int32_t)r;

    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)q/2;
}

void Decaps(uint32_t n, uint32_t q, const std::vector<int32_t>& t, const std::vector<int32_t>& s, const std::vector<int32_t>& c, std::vector<int32_t>& Hash_K, std::vector<int32_t> h_At, const std::vector<std::vector<int32_t>>& AT) {
    std::vector<int32_t> mprime(n, 0);
    std::vector<int32_t> u_j(n, 0);
    int32_t v_j = 0, dec_m = 0;

    for (uint32_t j = 0; j < n; ++j) {
        size_t off = (size_t)j * (n + 1);
        for (uint32_t i = 0; i < n; ++i) u_j[i] = c[off + i];
        v_j = c[off + n];
        Decrypt(v_j, u_j, s, q, dec_m);
        mprime[j] = (dec_m == (int32_t)(q/2)) ? 1 : 0;
    }

    std::vector<int32_t> K_cap = SHA256(concat(h_At, mprime));

    std::vector<int32_t> cchk; cchk.assign((size_t)n*(n+1), 0);
    std::vector<int32_t> u_tmp(n, 0);
    int32_t v_tmp = 0;

    for (uint32_t j = 0; j < n; ++j) {
        int32_t m_j_map = mprime[j] ? (int32_t)(q/2) : 0;
        std::vector<int32_t> r  = GenerateGaussianVector(n);
        std::vector<int32_t> e1 = GenerateGaussianVector(n);
        int32_t e2 = getRandomInt(-3, 3);

        Encrypt(n, q,
                const_cast<std::vector<int32_t>&>(t),
                u_tmp, v_tmp, (uint32_t)m_j_map, r, e1, e2, const_cast<std::vector<std::vector<int32_t>>&>(AT));

        size_t off = (size_t)j * (n + 1);
        for (uint32_t i = 0; i < n; ++i) cchk[off + i] = u_tmp[i];
        cchk[off + n] = v_tmp;
    }

    bool equal = (cchk.size() == c.size());
    if (equal) {
        for (size_t k = 0; k < c.size(); ++k)
            if (cchk[k] != c[k]) { equal = false; break; }
    }

    std::vector<int32_t> h_c = SHA256(c);
    if (equal) {
        Hash_K = SHA256( concat(K_cap, h_c) );
    } else {
        // implicit rejection 
    }
}

int main() {
    auto startTot = std::chrono::steady_clock::now();

    uint32_t n = 512;              
    uint32_t q = 3329;

    std::vector<std::vector<int32_t>> A;
    std::vector<int32_t> s, t;
    
    auto startK = std::chrono::steady_clock::now();
    KeyGen(n, q, A, s, t);
    auto endK = std::chrono::steady_clock::now();

    auto elapsedK = std::chrono::duration_cast<std::chrono::microseconds>(endK - startK);
    std::cout << "KeyGen time: " << elapsedK.count() << " mus\n";

    std::vector<std::vector<int32_t>> AT;
    transposeA(A, AT);

    std::vector<int32_t> m(n, 0);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> bit(0,1);
    for (uint32_t i = 0; i < n; ++i) m[i] = bit(rng);

    std::vector<int32_t> c;        
    std::vector<int32_t> K_enc; 
    
    std::vector<int32_t> At    = concat(flatten_matrix(A), t);
    std::vector<int32_t> h_At  = SHA256(At);

    auto startE = std::chrono::steady_clock::now();
    Encaps(n, q, t, m, c, K_enc, h_At, AT);
    auto endE = std::chrono::steady_clock::now();
    auto elapsedE = std::chrono::duration_cast<std::chrono::microseconds>(endE - startE);
    std::cout << "Encaps time: " << elapsedE.count() << " mus\n";


    std::vector<int32_t> z(32);
    for (int i = 0; i < 32; ++i) {
        z[i] = 42;
    }

    std::vector<int32_t> K_dec;

    auto startD = std::chrono::steady_clock::now();
    Decaps(n, q, t, s, c, K_dec, h_At, AT);
    auto endD = std::chrono::steady_clock::now();
    auto elapsedD = std::chrono::duration_cast<std::chrono::microseconds>(endD - startD);
    std::cout << "Decaps time: " << elapsedD.count() << " mus\n";
    
    bool sameK = (K_enc.size() == K_dec.size());
    if (sameK) {
        for (size_t i = 0; i < K_enc.size(); ++i) {
            if (K_enc[i] != K_dec[i]) { sameK = false; break; }
        }
    }
    
    auto endTot = std::chrono::steady_clock::now();

    auto elapsedTot = std::chrono::duration_cast<std::chrono::microseconds>(endTot - startTot);
    std::cout << "Tot time: " << elapsedTot.count() << " mus\n";
    return 0;
}
