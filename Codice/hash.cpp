#include <vector>
#include <cstdint>

#include <botan/hash.h>
#include <botan/hex.h>

#include "hash.h"

/////////////////////////////////////   SHAKE256    /////////////////////////////////////

/* 
    Non ho trovato un'implementazione ufficiale di SHA3, uso la libreria Botan che le implementa entrambe (forse c'è 
    qualcosa anche in OpenSSL)
*/

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
/* 
    SHAKE256 è una XOF:
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
    size_t total_byte  = 2 * n + 1;          // 513 totali

    // Byte totali necessari (per ogni bit mi serve total_bits, quindi 256 * total_bits)
    size_t out_bytes = msg_bits * total_byte;

    // Chiamo SHAKE-256 per ottenere out_bytes byte
    return shake256(seed_bytes, out_bytes);
    
}

// SHA3-256 usando Botan
std::vector<int32_t> SHA3_256(const std::vector<int32_t>& in)
{
    // converto il vettore di int32_t in byte (stesso schema di SHAKE256)
    std::vector<uint8_t> bytes;
    bytes.reserve(in.size() * 4);
    for (int32_t x : in)
    {
        uint32_t ux = static_cast<uint32_t>(x);
        bytes.push_back(static_cast<uint8_t>(ux & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 8)  & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 16) & 0xFF));
        bytes.push_back(static_cast<uint8_t>((ux >> 24) & 0xFF));
    }

    // creo l’oggetto hash SHA-3(256)
    auto hash = Botan::HashFunction::create_or_throw("SHA-3(256)");

    hash->update(bytes.data(), bytes.size());

    // alloco il buffer di output (32 byte) e lo riempio
    std::vector<uint8_t> digest(hash->output_length());
    hash->final(digest.data());

    // converto i byte in int32_t come facevi prima
    std::vector<int32_t> out(digest.size());
    for (size_t i = 0; i < digest.size(); ++i)
        out[i] = static_cast<int32_t>(digest[i]);

    return out; // 32 elementi → 256 bit
}
