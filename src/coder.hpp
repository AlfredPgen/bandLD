// Adaptive binary range coder (LZMA style) and an integer model for LD-band values.
//
// A signed integer v is coded as: bit length of |v| (0..63) through an adaptive 6-level bit tree,
// the sign, the two bits below the leading one with adaptive models, and the remaining bits raw.
// Every pair selects one of NCTX models from its context, so that values with a similar noise
// scale share statistics. All arithmetic is integer: encoder and decoder stay in step exactly.
#pragma once
#include <cstdint>
#include <cstring>
#include <vector>

namespace bcode {

constexpr int kProbBits = 11, kMoveBits = 5;
constexpr uint32_t kTop = 1u << 24;

struct Encoder {
    std::vector<uint8_t>& out;
    uint64_t low = 0; uint32_t range = 0xFFFFFFFFu; uint8_t cache = 0; uint64_t cache_size = 1;
    explicit Encoder(std::vector<uint8_t>& o) : out(o) {}
    inline void shift_low() {
        if ((uint32_t)low < 0xFF000000u || (low >> 32) != 0) {
            uint8_t c = cache;
            do { out.push_back(uint8_t(c + uint8_t(low >> 32))); c = 0xFF; } while (--cache_size != 0);
            cache = uint8_t(uint32_t(low) >> 24);
        }
        cache_size++;
        low = uint64_t(uint32_t(low) << 8);
    }
    inline void bit(uint16_t& p, uint32_t b) {
        const uint32_t bound = (range >> kProbBits) * p;
        if (!b) { range = bound; p = uint16_t(p + (((1u << kProbBits) - p) >> kMoveBits)); }
        else { low += bound; range -= bound; p = uint16_t(p - (p >> kMoveBits)); }
        while (range < kTop) { range <<= 8; shift_low(); }
    }
    inline void direct(uint64_t v, int nbits) {      // nbits raw bits of v, most significant first
        for (int i = nbits - 1; i >= 0; i--) {
            range >>= 1;
            if ((v >> i) & 1) low += range;
            while (range < kTop) { range <<= 8; shift_low(); }
        }
    }
    void finish() { for (int i = 0; i < 5; i++) shift_low(); }
};

struct Decoder {
    const uint8_t* p; const uint8_t* e;
    uint32_t range = 0xFFFFFFFFu, code = 0;
    Decoder(const uint8_t* b, size_t n) : p(b), e(b + n) {
        next();                                       // first byte is always 0
        for (int i = 0; i < 4; i++) code = (code << 8) | next();
    }
    inline uint8_t next() { return p < e ? *p++ : 0; }
    inline uint32_t bit(uint16_t& pr) {
        const uint32_t bound = (range >> kProbBits) * pr;
        uint32_t b;
        if (code < bound) { range = bound; pr = uint16_t(pr + (((1u << kProbBits) - pr) >> kMoveBits)); b = 0; }
        else { code -= bound; range -= bound; pr = uint16_t(pr - (pr >> kMoveBits)); b = 1; }
        while (range < kTop) { range <<= 8; code = (code << 8) | next(); }
        return b;
    }
    inline uint64_t direct(int nbits) {
        uint64_t v = 0;
        for (int i = 0; i < nbits; i++) {
            range >>= 1;
            uint32_t b = code >= range;
            if (b) code -= range;
            while (range < kTop) { range <<= 8; code = (code << 8) | next(); }
            v = (v << 1) | b;
        }
        return v;
    }
};

struct IntModel {
    uint16_t tree[64], sign, m1[64], m2[64][2];
    void reset() {
        const uint16_t h = 1u << (kProbBits - 1);
        for (auto& x : tree) x = h;
        sign = h;
        for (auto& x : m1) x = h;
        for (auto& x : m2) { x[0] = h; x[1] = h; }
    }
};

inline void put(Encoder& rc, IntModel& m, int64_t v) {
    const uint64_t u = v < 0 ? uint64_t(0) - uint64_t(v) : uint64_t(v);
    const int nb = u ? 64 - __builtin_clzll(u) : 0;
    uint32_t node = 1;
    for (int b = 5; b >= 0; b--) { const uint32_t t = (uint32_t(nb) >> b) & 1; rc.bit(m.tree[node], t); node = (node << 1) | t; }
    if (!nb) return;
    rc.bit(m.sign, v < 0);
    if (nb < 2) return;
    const uint32_t b1 = uint32_t(u >> (nb - 2)) & 1;
    rc.bit(m.m1[nb], b1);
    if (nb < 3) return;
    rc.bit(m.m2[nb][b1], uint32_t(u >> (nb - 3)) & 1);
    if (nb > 3) rc.direct(u & ((uint64_t(1) << (nb - 3)) - 1), nb - 3);
}

inline int64_t get(Decoder& rc, IntModel& m) {
    uint32_t node = 1;
    for (int b = 0; b < 6; b++) node = (node << 1) | rc.bit(m.tree[node]);
    const int nb = int(node - 64);
    if (!nb) return 0;
    const uint32_t neg = rc.bit(m.sign);
    uint64_t u = 1;
    if (nb >= 2) {
        const uint32_t b1 = rc.bit(m.m1[nb]);
        u = (u << 1) | b1;
        if (nb >= 3) {
            u = (u << 1) | rc.bit(m.m2[nb][b1]);
            if (nb > 3) u = (u << (nb - 3)) | rc.direct(nb - 3);
        }
    }
    return neg ? int64_t(uint64_t(0) - u) : int64_t(u);
}

// contexts: 256 noise-scale classes x 4 distance classes
constexpr int NSCALE = 256, NDIST = 4, NCTX = NSCALE * NDIST;
inline int dist_class(size_t d) { return d < 64 ? 0 : (d < 512 ? 1 : (d < 4096 ? 2 : 3)); }

struct Models {
    std::vector<IntModel> m;
    Models() : m(NCTX) {}
    void reset() { for (auto& x : m) x.reset(); }
    IntModel& at(int scale, size_t d) { return m[size_t(scale) * NDIST + dist_class(d)]; }
};

}  // namespace bcode
