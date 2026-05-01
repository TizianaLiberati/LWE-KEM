#pragma once
#include <cstdint>
#include <cstring>
#include <array>
#include <string>
#include <vector>

namespace simple_sha256 {

// --- util ---
static inline uint32_t rotr(uint32_t x, uint32_t n){ return (x >> n) | (x << (32U - n)); }
static inline uint32_t Ch(uint32_t x, uint32_t y, uint32_t z){ return (x & y) ^ (~x & z); }
static inline uint32_t Maj(uint32_t x, uint32_t y, uint32_t z){ return (x & y) ^ (x & z) ^ (y & z); }
static inline uint32_t BSIG0(uint32_t x){ return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22); }
static inline uint32_t BSIG1(uint32_t x){ return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25); }
static inline uint32_t SSIG0(uint32_t x){ return rotr(x,7) ^ rotr(x,18) ^ (x >> 3); }
static inline uint32_t SSIG1(uint32_t x){ return rotr(x,17) ^ rotr(x,19) ^ (x >> 10); }

class SHA256 {
public:
    SHA256(){ reset(); }

    void reset(){
        // FIPS 180-4, sezione 5.3.3
        h_[0]=0x6a09e667u; h_[1]=0xbb67ae85u; h_[2]=0x3c6ef372u; h_[3]=0xa54ff53au;
        h_[4]=0x510e527fu; h_[5]=0x9b05688cu; h_[6]=0x1f83d9abu; h_[7]=0x5be0cd19u;
        total_len_ = 0;
        buff_len_  = 0;
    }

    void update(const uint8_t* data, size_t len){
        total_len_ += len;
        // usa buffer interno per processare blocchi da 64 byte
        if (buff_len_ > 0){
            size_t to_copy = (len < (64 - buff_len_)) ? len : (64 - buff_len_);
            std::memcpy(buffer_ + buff_len_, data, to_copy);
            buff_len_ += to_copy;
            data += to_copy;
            len  -= to_copy;
            if (buff_len_ == 64){
                process_block(buffer_);
                buff_len_ = 0;
            }
        }
        while (len >= 64){
            process_block(data);
            data += 64;
            len  -= 64;
        }
        if (len > 0){
            std::memcpy(buffer_, data, len);
            buff_len_ = len;
        }
    }

    void update(const std::vector<uint8_t>& v){
        if(!v.empty()) update(v.data(), v.size());
    }
    void update(const std::string& s){
        if(!s.empty()) update(reinterpret_cast<const uint8_t*>(s.data()), s.size());
    }

    std::array<uint8_t,32> final(){
        // padding: 0x80, poi zeri fino a 56 mod 64, poi length (64 bit big-endian)
        uint8_t pad[64] = {0x80};
        size_t pad_len = (buff_len_ < 56) ? (56 - buff_len_) : (120 - buff_len_);
        update(pad, pad_len);

        uint64_t bit_len = static_cast<uint64_t>(total_len_) * 8ULL;
        uint8_t len_be[8];
        for(int i=0;i<8;i++) len_be[7-i] = static_cast<uint8_t>(bit_len >> (i*8));
        update(len_be, 8);

        std::array<uint8_t,32> out{};
        for(int i=0;i<8;i++){
            out[i*4+0] = static_cast<uint8_t>(h_[i] >> 24);
            out[i*4+1] = static_cast<uint8_t>(h_[i] >> 16);
            out[i*4+2] = static_cast<uint8_t>(h_[i] >> 8);
            out[i*4+3] = static_cast<uint8_t>(h_[i]);
        }
        // pronto per riuso
        reset();
        return out;
    }

    static std::array<uint8_t,32> hash(const uint8_t* data, size_t len){
        SHA256 ctx; ctx.update(data, len); return ctx.final();
    }
    static std::array<uint8_t,32> hash(const std::vector<uint8_t>& v){
        return hash(v.data(), v.size());
    }
    static std::array<uint8_t,32> hash(const std::string& s){
        return hash(reinterpret_cast<const uint8_t*>(s.data()), s.size());
    }

private:
    void process_block(const uint8_t* block){
        static const uint32_t K[64] = {
            0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
            0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
            0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
            0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
            0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
            0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
            0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
            0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
        };

        uint32_t W[64];
        // carica 16 parole big-endian
        for(int i=0;i<16;i++){
            W[i] = (static_cast<uint32_t>(block[i*4+0]) << 24)
                 | (static_cast<uint32_t>(block[i*4+1]) << 16)
                 | (static_cast<uint32_t>(block[i*4+2]) << 8)
                 | (static_cast<uint32_t>(block[i*4+3]));
        }
        // estendi a 64 parole
        for(int i=16;i<64;i++){
            W[i] = SSIG1(W[i-2]) + W[i-7] + SSIG0(W[i-15]) + W[i-16];
        }

        uint32_t a=h_[0], b=h_[1], c=h_[2], d=h_[3], e=h_[4], f=h_[5], g=h_[6], h=h_[7];

        for(int i=0;i<64;i++){
            uint32_t T1 = h + BSIG1(e) + Ch(e,f,g) + K[i] + W[i];
            uint32_t T2 = BSIG0(a) + Maj(a,b,c);
            h = g; g = f; f = e; e = d + T1;
            d = c; c = b; b = a; a = T1 + T2;
        }

        h_[0]+=a; h_[1]+=b; h_[2]+=c; h_[3]+=d; h_[4]+=e; h_[5]+=f; h_[6]+=g; h_[7]+=h;
    }

    uint32_t h_[8];
    uint8_t  buffer_[64];
    size_t   buff_len_ = 0;
    size_t   total_len_ = 0; // in byte
};

// --- API helper di alto livello ---

// digest in byte (32)
inline std::vector<uint8_t> sha256_bytes(const std::vector<uint8_t>& data){
    auto d = SHA256::hash(data);
    return std::vector<uint8_t>(d.begin(), d.end());
}
inline std::vector<uint8_t> sha256_bytes(const std::string& s){
    auto d = SHA256::hash(s);
    return std::vector<uint8_t>(d.begin(), d.end());
}

// digest in hex
inline std::string to_hex(const std::vector<uint8_t>& bytes){
    static const char* hex = "0123456789abcdef";
    std::string out; out.resize(bytes.size()*2);
    for(size_t i=0;i<bytes.size();++i){
        out[2*i  ] = hex[(bytes[i] >> 4) & 0xF];
        out[2*i+1] = hex[ bytes[i]       & 0xF];
    }
    return out;
}
inline std::string sha256_hex(const std::vector<uint8_t>& data){ return to_hex(sha256_bytes(data)); }
inline std::string sha256_hex(const std::string& s){ return to_hex(sha256_bytes(s)); }

// --- Adapter per il tuo codice che usa vector<int32_t> ---

// Serializza un vettore di int32_t in byte (little-endian 4 byte per elemento)
inline std::vector<uint8_t> i32_to_le_bytes(const std::vector<int32_t>& v){
    std::vector<uint8_t> out; out.reserve(v.size()*4);
    for(size_t i=0;i<v.size();++i){
        uint32_t x = static_cast<uint32_t>(v[i]);
        out.push_back((uint8_t)(x      & 0xFF));
        out.push_back((uint8_t)((x>>8) & 0xFF));
        out.push_back((uint8_t)((x>>16)& 0xFF));
        out.push_back((uint8_t)((x>>24)& 0xFF));
    }
    return out;
}

// Hash che accetta vector<int32_t> e restituisce digest come vector<int32_t> (32 elementi 0..255)
inline std::vector<int32_t> Hash_i32digest_from_i32input(const std::vector<int32_t>& in){
    auto bytes   = i32_to_le_bytes(in);
    auto digestB = sha256_bytes(bytes); // 32 byte
    std::vector<int32_t> out(32);
    for(size_t i=0;i<32;++i) out[i] = static_cast<int32_t>(digestB[i]);
    return out;
}

// Variante che restituisce direttamente i byte
inline std::vector<uint8_t> Hash_bytes_from_i32input(const std::vector<int32_t>& in){
    auto bytes = i32_to_le_bytes(in);
    return sha256_bytes(bytes);
}

} // namespace simple_sha256