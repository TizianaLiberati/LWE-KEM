#include <vector>
#include <array>
#include <cstdint>
#include <random>
#include <limits>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include "sha256.h"


/* Non ho trovato un'implementazione ufficiale di SHA3, uso la libreria Botan che le implementa entrambe (forse c'è qualcosa anche in OpenSSL)*/

#include <botan/hash.h>
#include <botan/hex.h>

#include <memory>

using simple_sha256::i32_to_le_bytes;
using simple_sha256::sha256_bytes;


/////////////////////////////////////   SHAKE256    /////////////////////////////////////

// Funzioni ausiliarie per SHAKE256

// Botan vuole byte (uint8_t), noi lavoriamo in int32_t
std::vector<uint8_t> int32_to_bytes(const std::vector<int32_t>& v)
{
    std::vector<uint8_t> out;
    out.reserve(v.size());
    for (int32_t x : v) {
        out.push_back(static_cast<uint8_t>(x & 0xFF));
    }
    return out;
}

// Prende in input i byte del seed e in output restituisce il numero di byte necessari 
std::vector<uint8_t> shake256(const std::vector<uint8_t>& input, size_t out_len_bytes)
{
    const size_t out_bits = out_len_bytes * 8; // Usa i bit solo nel nome (Es: "SHAKE-256(520)" se out_len_bytes = 65)

    const std::string xof = "SHAKE-256(" + std::to_string(out_bits) + ")"; // uso l'algoritmo SHAKE256

    auto shake = Botan::HashFunction::create_or_throw(xof);

    shake->update(input.data(), input.size()); // inserisco il seed e la sua lunghezza
    std::vector<uint8_t> out(out_len_bytes); // creo vettore che riempirò

    shake->final(out.data());

    return out;
}

// Prende in input il seed (int32_t), lo converte in byte (uint8_t) e in output restituisce la stringa dei coins da cui ricaverò i noise
/* SHAKE256 è una XOF:
    - stesso input -> stesso output
    - scelta lunghezza dell'output
*/
std::vector<uint8_t> xof_coins(const std::vector<int32_t>& seed_int32, size_t n, size_t msg_bits)
{
    // const size_t n = 512;
    // const size_t msg_bits = 256;

    // Converto il seed (bits) in byte
    std::vector<uint8_t> seed_bytes = int32_to_bytes(seed_int32);

    // Calcolo quanti byte servono 
    // Ogni bit del messaggio usa:
    // - 256 byte per r
    // - 256 byte per e1
    // - 1 byte per e2
    size_t total_bits  = 2 * n + 1;          // 513 totali

    // Byte totali necessari (per ogni bit mi serve total_bits, quindi 256 * total_bits)
    size_t out_bytes = msg_bits * total_bits;

    // Chiamo SHAKE-256 per ottenere out_bytes byte
    return shake256(seed_bytes, out_bytes);
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int mod(int a, int b)
{
    return (a % b + b) % b;
}

// TODO: sostituire con SHA3-256
static std::vector<int32_t> SHA256(const std::vector<int32_t> &in)
{
    std::vector<uint8_t> bytes;
    bytes.reserve(in.size() * 4);
    for (int32_t x : in)
    {
        uint32_t ux = (uint32_t)x;
        bytes.push_back((uint8_t)(ux & 0xFF));
        bytes.push_back((uint8_t)((ux >> 8) & 0xFF));
        bytes.push_back((uint8_t)((ux >> 16) & 0xFF));
        bytes.push_back((uint8_t)((ux >> 24) & 0xFF));
    }
    auto digest = simple_sha256::sha256_bytes(bytes); // 32 byte
    std::vector<int32_t> out(32);
    for (size_t i = 0; i < 32; ++i)
        out[i] = (int32_t)digest[i];
    return out;
}

/////////////////////////////////////   Funzioni ausiliarie /////////////////////////////////////

std::vector<int32_t> concat(const std::vector<int32_t> &a, const std::vector<int32_t> &b)
{
    std::vector<int32_t> out;
    out.reserve(a.size() + b.size());
    for (size_t i = 0; i < a.size(); ++i)
        out.push_back(a[i]);
    for (size_t i = 0; i < b.size(); ++i)
        out.push_back(b[i]);
    return out;
}

std::vector<int32_t> flatten_matrix(const std::vector<std::vector<int32_t>> &A)
{
    std::vector<int32_t> out;
    if (A.empty())
        return out;
    const size_t n = A.size();
    out.reserve(n * A[0].size());
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            out.push_back(A[i][j]);
    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Generazioni in KeyGen   /////////////////////////////////////

// sample_discrete_gaussian e GenerateGaussianVector le uso per la generazione dell'errore e nella funzione KeyGen
int32_t sample_discrete_gaussian(double sigma)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::normal_distribution<double> norm(0.0, sigma);

    double x = norm(gen);

    double bound = 6.0 * sigma;
    if (x > bound)
        x = bound;
    if (x < -bound)
        x = -bound;

    long long k = std::llround(x);

    if (k > std::numeric_limits<int32_t>::max())
        k = std::numeric_limits<int32_t>::max();
    if (k < std::numeric_limits<int32_t>::min())
        k = std::numeric_limits<int32_t>::min();
    return static_cast<int32_t>(k);
}

std::vector<int32_t> GenerateGaussianVector(size_t n)
{
    std::vector<int32_t> vec;
    vec.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        vec.push_back(sample_discrete_gaussian(2.3));
    }
    return vec;
}

// getRandomInt la useremo solo per campionare z
int32_t getRandomInt(int min, int max)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(min, max);
    return distr(gen);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Generazione noise r, e1, e2 /////////////////////////////////////

// sample_e2 restituisce un vettore in [-3, 3] (binomiale centrata)
/*
- coins: output di SHAKE256
- pos: posizione dei coins presi della stringa XOF
*/
int32_t sample_e2(uint8_t coins)
{
    uint32_t v = static_cast<uint32_t>(coins) % 7;   // v \in [0, 6]
    return static_cast<int32_t>(v) - 3; // centriamo v \in [-3, 3]
}


/* Road Map:
    SHAKE -> Valori campionati da una gaussiana con deviazione standard = 2.3
    -> SHAKE restituisce un output che si comporta come un'uniforme DISCRETA tra [0, 1] ({0, 1})
    -> da uniforme discreta dobbiamo passare ad un'uniforme CONTINUA sempre in [0, 1] (u \in (0, 1)) quindi normalizziamo
    -> CDF inversa prende in input un valore reale campionato da un'uniforme CONTINUA (u \in (0, 1)) e restituisce un valore appartenente ad una normale (Z \ondina N(0, 1))

*/

// esistono altre tecniche per mappare un'uniforme in una gaussiana
// (es. Box Muller, Ziggurat, Ratio-of-uniforms) online ci sono pareri moolto diversi quindi si può pensare di usare
// un metodo diverso da CDF inversa

/*CDF inversa:
Il metodo si basa sul fatto che se X è una variabile casuale continua con una funzione di ripartizione strettamente 
crescente FX e Y = FX(X), allora Y ha una distribuzione uniforme nell'intervallo [FX_min , FX_max].

Il metodo procede come segue:

data una variabile casuale uniforme continua U in [0, 1] e una funzione di ripartizione invertibile F, 
la variabile casuale X = F^{−1}(U) è distribuita secondo F (o, equivalentemente X ha la distribuzione F).
*/

// Normalizzazione: dai valori della SHAKE (b) resituiamo valori interi (double) da un'uniforme continua (u)
double uniform01_from_byte(uint8_t b)
{
    // b \in [0, 255] -> u = (b+0.5)/256 \in [0,1]
    double u = (static_cast<double>(b) + 0.5) / 256.0; 

    // il controllo serve per non avere u = 0 o u = 1 perchè la CDF tende a inf (+ o -) quando u = 1 o 0
    /*
    const double eps = 1e-12;
    if (u <= 0.0) u = eps; // in teoria non dovremmo avere queste casistiche
    if (u >= 1.0) u = 1.0 - eps;
    */

    return u;
}

// implementa CDF della N(0,1) (gaussiana standard)
double normal_cdf(double x)
{
    // funzione per calcolar la CDF (da teoria)
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Implementa CDF inversa (gaussiana standard), approssimiamo con ricerca binaria
/* 
   In teoria si usano formule approssimate ma (da quello che ho visto in giro) si semplifica cercando x in 
   un intervallo tale che CDF(x) sia vicino a u.
*/
double normal_icdf(double u)
{
    // facciamo comunque il controllo che u sia diverso da 0 e 1
    const double eps = 1e-12;
    if (u <= 0.0) u = eps;
    if (u >= 1.0) u = 1.0 - eps;

    // poniamo come intervallo di ricerca per x, t.c CDF(x) = u, [-10, 10]
    // per una gaussiana standard la probabilità P(|Z| > 10) è dell'ordine di 10^{-23} (tabelle)
    double lo = -10.0;
    double hi =  10.0;

    for (int iter = 0; iter < 60; ++iter)
    {
        double mid = 0.5 * (lo + hi); // all'inizio mid = 0 poi itero
        double c = normal_cdf(mid); // c = CDF(mid)

        if (c < u)
            lo = mid;  // serve un x più grande, valori più grandi
        else
            hi = mid;  // serve un x più piccolo, valori più piccoli
    }

    // abbiamo iterato quindi lo, hi diversi da 0
    return 0.5 * (lo + hi);
}

// SHAKE -> N(0,1) come double (unisce le cose fatte prima)
double sample_normal01_from_byte(uint8_t coins)
{
    double u = uniform01_from_byte(coins); // uniforme u(0,1)
    double z = normal_icdf(u);                 // Z ~ N(0,1) tramite la CDF inversa
    return z;
}

// Funzione "finale": restituisce valori campionati da una gaussiana di deviazione standard \sigma = 2.3
int32_t sample_discrete_gaussian_sigma_from_byte(uint8_t coins, double sigma)
{
    // 1) Z ~ N(0,1) preso dai byte della SHAKE
    double z = sample_normal01_from_byte(coins);

    // 2) scala con sigma: ora val ~ N(0, sigma^2)
    double val = sigma * z;

    // 3) tronchiamo a ±6\sigma (evita valori troppo grandi, in una gaussiana i valori sulle code sono rari, 
    // aka le code sono "basse")
    double bound = 6.0 * sigma;
    if (val >  bound) val =  bound;
    if (val < -bound) val = -bound;

    // 4) arrotonda all'intero più vicino
    long long k = std::llround(val);

    // 5) controllo caso di overflow (non dovrebbe mai scattare in pratica)
    // avendo troncato a +- 6 \sigma dovrei avere arrotondando al massimo k = +- 14
    // in teoria non dovrei uscire da int32_t ma prima di commentare faccio qualche prova 
    if (k > std::numeric_limits<int32_t>::max())
        k = std::numeric_limits<int32_t>::max();
    if (k < std::numeric_limits<int32_t>::min())
        k = std::numeric_limits<int32_t>::min();

    // 6) restituisci come int32_t
    return static_cast<int32_t>(k);
}

// Creiamo una struttura per i noise che useremo spesso
struct NoiseTriple
{
    std::vector<int32_t> r;
    std::vector<int32_t> e1;
    int32_t e2;
};

// Funzione che genera i noise per un bit del messaggio
NoiseTriple GenerateNoisesForOneBit(const std::vector<uint8_t>& coins, size_t& pos, const size_t n)
{
    // const size_t n = 512;
    double sigma = 2.3; 

    NoiseTriple nt;

    // alloco lunghezza noise r ed e1
    nt.r.reserve(n); 
    nt.e1.reserve(n);

    // r: 256 campioni gaussiani
    for (size_t j = 0; j < n; ++j)
    {
        uint8_t b = coins[pos];
        ++pos;
        nt.r.push_back(sample_discrete_gaussian_sigma_from_byte(b, sigma));
    }

    // e1: 256 campioni gaussiani
    for (size_t j = 0; j < n; ++j)
    {
        uint8_t b = coins[pos];
        ++pos;
        nt.e1.push_back(sample_discrete_gaussian_sigma_from_byte(b, sigma));
    }

    // e2: 1 byte -> [-3,3]
    {
        uint8_t b = coins[pos];
        ++pos;
        nt.e2 = sample_e2(b);
    }

    return nt;   
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Generazione matrice A e calcolo trasposta di A   /////////////////////////////////////
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

void transposeA(const std::vector<std::vector<int32_t>> &A, std::vector<std::vector<int32_t>> &AT)
{
    const size_t n = A.size();
    AT.assign(n, std::vector<int32_t>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            AT[j][i] = A[i][j];
}

/////////////////////////////////////   Calcolo segreto s in KeyGen /////////////////////////////////////

int32_t sample_eta_centered_binomial(uint8_t eta, std::mt19937 &gen)
{
    std::uniform_int_distribution<uint8_t> dis(0, 1);
    uint8_t sum1 = 0, sum2 = 0;
    for (uint8_t i = 0; i < eta; ++i)
    {
        sum1 += dis(gen);
        sum2 += dis(gen);
    }
    return (int32_t)sum1 - (int32_t)sum2;
}

std::vector<int32_t> sample_vector_binomial(uint32_t n)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<int32_t> result;
    result.reserve(n);
    for (uint32_t i = 0; i < n; ++i)
        result.push_back(sample_eta_centered_binomial(3, gen));
    return result;
}

/*
    https://github.com/pq-crystals/dilithium/blob/master/ref/randombytes.c
    qui generano bit random in modo crittograficamente sicuro, nel mio caso mt19937 è random ma non compatibile con richieste di sicurezza in crittografia
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   KeyGen  /////////////////////////////////////

void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &s_k, std::vector<int32_t> &t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, q);

    std::vector<int32_t> z(256); //256 bits
    for (int i = 0; i < 256; ++i)
        z[i] = getRandomInt(0, 1);

    // secret key s_k data dalla concat di s e z
    s_k = concat(s, z);

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Encrypt /////////////////////////////////////
void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2, const std::vector<std::vector<int32_t>> &AT)
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Encaps  /////////////////////////////////////

void Encaps(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &c, const std::vector<std::vector<int32_t>> &A, const std::vector<std::vector<int32_t>> &AT, std::vector<int32_t> &Hash_K)
{
    const size_t msg_bits = 256;

    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh = SHA256(At); // potrei far in modo che SHA256 restituisca in uint8_t (magari se passiamo in SHA3-526)

    // Genero il plaintext lungo 256 random
    std::vector<int32_t> m(msg_bits, 0);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> bit(0, 1);
    for (uint32_t i = 0; i < msg_bits; ++i)
        m[i] = bit(rng);

    // Creo il seed che userò nella XOF
    std::vector<int32_t> pkh_m = concat(pkh, m);
    std::vector<int32_t> seed = SHA256(pkh_m);


    // implementa XOF e prendi i seed per i noise
    std::vector<uint8_t> coins = xof_coins(seed, n, msg_bits);

    // std::vector<int32_t> r = GenerateGaussianVector(n);
    // std::vector<int32_t> e1 = GenerateGaussianVector(n);

    // int32_t e2 = getRandomInt((-1) * 3, 3); 

    // TODO: RIVEDI COME GENERARE K_cap
    std::vector<int32_t> K_cap = SHA256(concat(pkh, m));

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

    std::vector<int32_t> h_c = SHA256(c);
    Hash_K = SHA256(concat(K_cap, h_c));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Decrypt /////////////////////////////////////

void Decrypt(int32_t v_i, const std::vector<int32_t> &u, const std::vector<int32_t> &s_k, uint32_t q, int32_t &decrypt_i)
{
    long long dot = 0;

    // Considero solo la prima parte di s_k, quella relativa ad s
    const size_t n = u.size(); // la dimensione di s (senza z) è n = dimensione di u
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n); // definisco s uguale al primo elemento di s_k fino all'elemento n-esimo di s_k

    for (size_t i = 0; i < s.size(); ++i)
        dot += (long long)s[i] * (long long)u[i];

    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0)
        r += (long long)q;
    int32_t mu = (int32_t)r;

    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)q / 2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Decaps  /////////////////////////////////////

void Decaps(uint32_t n, uint32_t q, const std::vector<int32_t> &t, const std::vector<int32_t> &s_k, const std::vector<int32_t> &c, std::vector<int32_t> &Hash_K, const std::vector<std::vector<int32_t>> &A, const std::vector<std::vector<int32_t>> &AT)
{
    const size_t msg_bits = 256;

    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh = SHA256(At);

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

    std::vector<int32_t> K_cap = SHA256(concat(pkh, mprime));

    std::vector<int32_t> cchk;
    cchk.assign((size_t)msg_bits * (n + 1), 0);
    std::vector<int32_t> u_tmp(n, 0);
    int32_t v_tmp = 0;

    // Creo il seed che userò nella XOF
    std::vector<int32_t> pkh_mprime = concat(pkh, mprime);
    std::vector<int32_t> seed = SHA256(pkh_mprime);


    // implementa XOF e prendi i seed per i noise
    std::vector<uint8_t> coins = xof_coins(seed, n, msg_bits);
    size_t pos = 0; // indice dentro coins

    for (uint32_t j = 0; j < msg_bits; ++j)
    {
        int32_t m_j_map = mprime[j] ? (int32_t)(q / 2) : 0;
        // std::vector<int32_t> r = GenerateGaussianVector(n);
        // std::vector<int32_t> e1 = GenerateGaussianVector(n);
        // int32_t e2 = getRandomInt(-3, 3);

        // genera r, e1, e2 per ogni bit del plaintext (nel for)
        NoiseTriple noise_i = GenerateNoisesForOneBit(coins, pos, n);

        std::vector<int32_t>& r  = noise_i.r;
        std::vector<int32_t>& e1 = noise_i.e1;
        int32_t e2 = noise_i.e2;
        Encrypt(n, q, const_cast<std::vector<int32_t> &>(t), u_tmp, v_tmp, (uint32_t)m_j_map, r, e1, e2, const_cast<std::vector<std::vector<int32_t>> &>(AT));

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
    
    std::vector<int32_t> h_c = SHA256(c);
    if (equal)
    {
        Hash_K = SHA256(concat(K_cap, h_c));
        // std::cout << "ha funzionato\n";
    }
    else
    {
        // implicit rejection
        std::vector<int32_t> z(s_k.begin() + n, s_k.end());
        Hash_K = SHA256(concat(z, h_c));
        // std::cout << "non ha funzionato\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    std::vector<int32_t> K_enc;

    /* Lo sposto in Encaps

    std::vector<int32_t> At = concat(flatten_matrix(A), t);
    std::vector<int32_t> pkh_At = SHA256(At);
    */
  
    auto startE = std::chrono::steady_clock::now();
    Encaps(n, q, t, c, A, AT, K_enc); // K_enc ha 32 byte = 256 bit
    auto endE = std::chrono::steady_clock::now();
    auto elapsedE = std::chrono::duration_cast<std::chrono::microseconds>(endE - startE);
    std::cout << "Encaps time: " << elapsedE.count() << " mus\n";

    std::vector<int32_t> K_dec;

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