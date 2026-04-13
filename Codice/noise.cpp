#include <cstddef>
#include <cmath>     // std::sqrt, std::erfc, std::llround
#include <limits>
#include <vector>
#include <cstdint>

#include "config.h"
#include "noise.h"

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

/*
CDF inversa:
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

// Funzione che genera i noise per un bit del messaggio
NoiseTriple GenerateNoisesForOneBit(const std::vector<uint8_t>& coins,
                                    size_t& pos,
                                    const size_t n,
                                    const LweKemConfig& cfg)
{
    // const size_t n = 512;
    // double sigma = 2.3; 

    NoiseTriple nt;

    // alloco lunghezza noise r ed e1
    nt.r.reserve(n); 
    nt.e1.reserve(n);

    // r: 256 campioni gaussiani
    for (size_t j = 0; j < n; ++j)
    {
        uint8_t b = coins[pos];
        ++pos;
        nt.r.push_back(sample_from_kind(b, cfg.r_kind, cfg));
    }

    // e1: 256 campioni gaussiani
    for (size_t j = 0; j < n; ++j)
    {
        uint8_t b = coins[pos];
        ++pos;
        nt.e1.push_back(sample_from_kind(b, cfg.e1_kind, cfg));
    }

    // e2: 1 byte -> [-3,3]
    {
        uint8_t b = coins[pos];
        ++pos;
        nt.e2 = sample_e2_from_kind(b, cfg);
    }

    return nt;   
}

int32_t sample_from_kind(uint8_t coins, NoiseKind kind, const LweKemConfig& cfg)
{
    switch (kind) {
        case NoiseKind::BINOMIAL_ETA1:
            return static_cast<int32_t>(coins & 1) - static_cast<int32_t>((coins >> 1) & 1);

        case NoiseKind::BINOMIAL_ETA3: {
            int32_t a = static_cast<int32_t>(coins & 1)
                      + static_cast<int32_t>((coins >> 1) & 1)
                      + static_cast<int32_t>((coins >> 2) & 1);
            int32_t b = static_cast<int32_t>((coins >> 3) & 1)
                      + static_cast<int32_t>((coins >> 4) & 1)
                      + static_cast<int32_t>((coins >> 5) & 1);
            return a - b;
        }

        case NoiseKind::GAUSSIAN:
            return sample_discrete_gaussian_sigma_from_byte(coins, cfg.gauss_sigma);

        case NoiseKind::UNIFORM_CENTERED_3:
            return sample_e2(coins);
    }

    return 0;
}

int32_t sample_e2_from_kind(uint8_t coins, const LweKemConfig& cfg)
{
    return sample_from_kind(coins, cfg.e2_kind, cfg);
}