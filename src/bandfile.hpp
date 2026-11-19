// bandLD file format (.bld), one file per chromosome.
//
// SNPs are sorted by position. Row j holds one value for every earlier SNP i with
// pos[j] - pos[i] <= window_bp, i ascending ("backward" rows: a row is complete as soon as SNP j
// has been processed, and a region [a, b) is the rows a..b-1 restricted to columns >= a).
//
// Value of pair (i, j):
//   exact mode (round_bits = 0):  D = C' - floor((sx'_i * sx'_j + N/2) / N), where C' = sum_s k'_i[s] k'_j[s]
//       is the exact integer cross product of the reduced dosage numerators k' = k / g (g = the
//       SNP's common divisor of all its numerators: 255 for 8-bit hard calls, usually 1 for
//       dosages) and sx' = sx / g. With the per-SNP sums this gives back r exactly as a double,
//       evaluated in the original units like fastLD:
//           C = C' g_i g_j,   r = (N*C - sx_i*sx_j) / sqrt(V_i * V_j),   V = N*sxx - sx^2
//   rounded mode (round_bits = B > 0):  q = round(r * 2^B)
// Values are entropy coded (coder.hpp) in independent chunks of rows_per_chunk rows.
//
// Layout:  header (HDR_BYTES) | row chunks | SNP table (zstd) | chunk index (u64 offsets, n_chunks+1)
#pragma once
#include "common.hpp"
#include "coder.hpp"

constexpr size_t HDR_BYTES = 160;
static const char BLD_MAGIC[8] = {'B', 'A', 'N', 'D', 'L', 'D', 0, 1};

struct BandHeader {
    uint32_t version = 1, round_bits = 0;
    uint64_t n_samples = 0, n_snps = 0, window_bp = 0;
    uint32_t rows_per_chunk = 64, reserved = 0;
    uint64_t n_chunks = 0, n_pairs = 0;
    uint64_t table_off = 0, table_csize = 0, table_usize = 0, index_off = 0;
    char chrom[32] = {0};
};

struct BandSnps {                     // SNP table, all vectors of length M
    vector<uint32_t> pos;
    vector<int64_t> sx, sxx;          // sums of the numerators and of their squares (original units)
    vector<uint32_t> gcd;             // common divisor of the SNP's numerators
    vector<uint8_t> bits;             // probability bits of the source record (numerators are in units of 1/(2^bits-1))
    vector<string> rsid, a1, a2;      // dosage counted: allele2
    size_t size() const { return pos.size(); }
    // derived
    vector<uint8_t> vbits;            // bit length of V / g^2 (noise-scale context, reduced units)
    vector<double> dV;                // double(V), original units
    int64_t sxr(size_t j) const { return sx[j] / int64_t(gcd[j]); }   // reduced sum
    vector<uint64_t> lo;              // first index i with pos[i] + window >= pos[j]
    void derive(uint64_t N, uint64_t window_bp);
    double maf(size_t j, uint64_t N) const {
        const double den = double((uint32_t(1) << bits[j]) - 1), f = double(sx[j]) / (2.0 * den * double(N));
        return f < 0.5 ? f : 1.0 - f;
    }
};

inline int bitlen128(unsigned __int128 v) {
    const uint64_t hi = uint64_t(v >> 64), lo = uint64_t(v);
    return hi ? 128 - __builtin_clzll(hi) : (lo ? 64 - __builtin_clzll(lo) : 0);
}
inline unsigned __int128 snp_V(int64_t sx, int64_t sxx, uint64_t N) {
    return (unsigned __int128)N * (unsigned __int128)sxx - (unsigned __int128)sx * (unsigned __int128)sx;
}
// expected cross product under independence, rounded: floor((sx_i * sx_j + N/2) / N)
inline int64_t expect_C(int64_t sxi, int64_t sxj, uint64_t N) {
    return int64_t(((unsigned __int128)sxi * (unsigned __int128)sxj + N / 2) / N);
}
inline int scale_ctx(const BandSnps& s, size_t i, size_t j) { return std::min(255, int(s.vbits[i]) + int(s.vbits[j])); }
// r from the exact reduced cross product (same expression as fastLD, in the original units)
inline double r_exact(int64_t Cr, const BandSnps& s, size_t i, size_t j, uint64_t N) {
    const __int128 C = (__int128)Cr * s.gcd[i] * s.gcd[j];
    const __int128 num = (__int128)N * C - (__int128)s.sx[i] * s.sx[j];
    return double(num) / sqrt(s.dV[i] * s.dV[j]);
}

// ---- file access
void write_header(FILE* f, const BandHeader& h);
struct BandFile {
    string path; int fd = -1; uint64_t fsize = 0;
    BandHeader h; BandSnps snps; vector<uint64_t> chunk_off;
    ~BandFile();
    void open(const string& fn);
    // decode chunk c; rows[j - first_row] = values for i = lo[j] .. j-1 (exact mode: reduced C', rounded mode: q)
    void read_chunk(size_t c, vector<vector<int64_t>>& rows) const;
};
vector<uint8_t> pack_snp_table(const BandSnps& s);        // uncompressed image
void encode_chunk(const BandHeader& h, const BandSnps& s, size_t first_row, const vector<vector<int64_t>>& rows,
                  size_t row0_in_rows, size_t nrows, vector<uint8_t>& out, bcode::Models& models);
