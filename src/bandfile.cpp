#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

#include "bandfile.hpp"

void BandSnps::derive(uint64_t N, uint64_t window_bp) {
    const size_t M = size();
    vbits.resize(M); dV.resize(M); lo.resize(M);
    size_t l = 0;
    for (size_t j = 0; j < M; j++) {
        const unsigned __int128 V = snp_V(sx[j], sxx[j], N);
        vbits[j] = (uint8_t)bitlen128(V / ((unsigned __int128)gcd[j] * gcd[j])); dV[j] = double(V);
        while (uint64_t(pos[l]) + window_bp < uint64_t(pos[j])) l++;
        lo[j] = l;
    }
}

template <class T> static void put_le(uint8_t*& p, T v) { memcpy(p, &v, sizeof(T)); p += sizeof(T); }
template <class T> static T get_le(const uint8_t*& p) { T v; memcpy(&v, p, sizeof(T)); p += sizeof(T); return v; }

void write_header(FILE* f, const BandHeader& h) {
    uint8_t buf[HDR_BYTES]; memset(buf, 0, sizeof buf);
    uint8_t* p = buf;
    memcpy(p, BLD_MAGIC, 8); p += 8;
    put_le(p, h.version); put_le(p, h.round_bits);
    put_le(p, h.n_samples); put_le(p, h.n_snps); put_le(p, h.window_bp);
    put_le(p, h.rows_per_chunk); put_le(p, h.reserved);
    put_le(p, h.n_chunks); put_le(p, h.n_pairs);
    put_le(p, h.table_off); put_le(p, h.table_csize); put_le(p, h.table_usize); put_le(p, h.index_off);
    memcpy(p, h.chrom, 32);
    if (fseeko(f, 0, SEEK_SET) != 0 || fwrite(buf, 1, HDR_BYTES, f) != HDR_BYTES) die("cannot write the band file header");
}

vector<uint8_t> pack_snp_table(const BandSnps& s) {
    const size_t M = s.size();
    size_t strbytes = 0;
    for (size_t j = 0; j < M; j++) strbytes += s.rsid[j].size() + s.a1[j].size() + s.a2[j].size() + 3;
    vector<uint8_t> out(M * (4 + 8 + 8 + 4 + 1) + strbytes);
    uint8_t* p = out.data();
    memcpy(p, s.pos.data(), M * 4); p += M * 4;
    memcpy(p, s.sx.data(), M * 8); p += M * 8;
    memcpy(p, s.sxx.data(), M * 8); p += M * 8;
    memcpy(p, s.gcd.data(), M * 4); p += M * 4;
    memcpy(p, s.bits.data(), M); p += M;
    for (size_t j = 0; j < M; j++)
        for (const string* t : {&s.rsid[j], &s.a1[j], &s.a2[j]}) { memcpy(p, t->c_str(), t->size() + 1); p += t->size() + 1; }
    return out;
}

BandFile::~BandFile() { if (fd >= 0) close(fd); }

static void pread_all(int fd, void* buf, size_t n, uint64_t off, const string& fn) {
    uint8_t* p = (uint8_t*)buf;
    while (n) {
        ssize_t r = pread(fd, p, n, (off_t)off);
        if (r <= 0) die("short read in '" + fn + "' (truncated file?)");
        p += r; n -= r; off += r;
    }
}

void BandFile::open(const string& fn) {
    path = fn;
    fd = ::open(fn.c_str(), O_RDONLY);
    if (fd < 0) die("cannot open band file '" + fn + "'");
    struct stat sb;
    if (fstat(fd, &sb) != 0) die("cannot stat '" + fn + "'");
    fsize = (uint64_t)sb.st_size;
    if (fsize < HDR_BYTES) die("'" + fn + "' is not a bandLD file (too short)");
    uint8_t buf[HDR_BYTES];
    pread_all(fd, buf, HDR_BYTES, 0, fn);
    if (memcmp(buf, BLD_MAGIC, 8) != 0) die("'" + fn + "' is not a bandLD file");
    const uint8_t* p = buf + 8;
    h.version = get_le<uint32_t>(p); h.round_bits = get_le<uint32_t>(p);
    h.n_samples = get_le<uint64_t>(p); h.n_snps = get_le<uint64_t>(p); h.window_bp = get_le<uint64_t>(p);
    h.rows_per_chunk = get_le<uint32_t>(p); h.reserved = get_le<uint32_t>(p);
    h.n_chunks = get_le<uint64_t>(p); h.n_pairs = get_le<uint64_t>(p);
    h.table_off = get_le<uint64_t>(p); h.table_csize = get_le<uint64_t>(p); h.table_usize = get_le<uint64_t>(p); h.index_off = get_le<uint64_t>(p);
    memcpy(h.chrom, p, 32); h.chrom[31] = 0;
    if (h.version != 1) die("'" + fn + "': unsupported format version " + std::to_string(h.version));
    if (h.table_off == 0) die("'" + fn + "' is incomplete (the build did not finish)");
    const size_t M = h.n_snps;
    if (h.table_off + h.table_csize > fsize || h.index_off + (h.n_chunks + 1) * 8 > fsize || h.rows_per_chunk == 0 ||
        h.n_chunks != (M + h.rows_per_chunk - 1) / h.rows_per_chunk || h.table_usize < M * 25)
        die("'" + fn + "' is corrupt (header offsets)");
    vector<uint8_t> ct(h.table_csize), ut(h.table_usize);
    pread_all(fd, ct.data(), ct.size(), h.table_off, fn);
    const size_t got = ZSTD_decompress(ut.data(), ut.size(), ct.data(), ct.size());
    if (ZSTD_isError(got) || got != ut.size()) die("'" + fn + "' is corrupt (SNP table)");
    const uint8_t* q = ut.data(); const uint8_t* qe = q + ut.size();
    snps.pos.resize(M); snps.sx.resize(M); snps.sxx.resize(M); snps.gcd.resize(M); snps.bits.resize(M);
    memcpy(snps.pos.data(), q, M * 4); q += M * 4;
    memcpy(snps.sx.data(), q, M * 8); q += M * 8;
    memcpy(snps.sxx.data(), q, M * 8); q += M * 8;
    memcpy(snps.gcd.data(), q, M * 4); q += M * 4;
    memcpy(snps.bits.data(), q, M); q += M;
    snps.rsid.resize(M); snps.a1.resize(M); snps.a2.resize(M);
    for (size_t j = 0; j < M; j++)
        for (string* t : {&snps.rsid[j], &snps.a1[j], &snps.a2[j]}) {
            const uint8_t* z = (const uint8_t*)memchr(q, 0, qe - q);
            if (!z) die("'" + fn + "' is corrupt (SNP strings)");
            t->assign((const char*)q, z - q); q = z + 1;
        }
    for (size_t j = 1; j < M; j++) if (snps.pos[j] < snps.pos[j - 1]) die("'" + fn + "' is corrupt (positions not sorted)");
    for (size_t j = 0; j < M; j++) if (snps.gcd[j] == 0 || snps.sx[j] % int64_t(snps.gcd[j])) die("'" + fn + "' is corrupt (SNP divisors)");
    snps.derive(h.n_samples, h.window_bp);
    chunk_off.resize(h.n_chunks + 1);
    if (h.n_chunks + 1) pread_all(fd, chunk_off.data(), (h.n_chunks + 1) * 8, h.index_off, fn);
    for (size_t c = 0; c < h.n_chunks; c++)
        if (chunk_off[c] > chunk_off[c + 1] || chunk_off[c + 1] > h.table_off) die("'" + fn + "' is corrupt (chunk index)");
}

void encode_chunk(const BandHeader& h, const BandSnps& s, size_t first_row, const vector<vector<int64_t>>& rows,
                  size_t row0_in_rows, size_t nrows, vector<uint8_t>& out, bcode::Models& models) {
    out.clear();
    models.reset();
    bcode::Encoder rc(out);
    const bool exact = h.round_bits == 0;
    const uint64_t N = h.n_samples;
    for (size_t t = 0; t < nrows; t++) {
        const size_t j = first_row + t, l = s.lo[j];
        const vector<int64_t>& row = rows[row0_in_rows + t];
        for (size_t i = l; i < j; i++) {
            const int64_t v = row[i - l];
            if (exact) bcode::put(rc, models.at(scale_ctx(s, i, j), j - i), v - expect_C(s.sxr(i), s.sxr(j), N));
            else bcode::put(rc, models.at(scale_ctx(s, i, j), j - i), v);
        }
    }
    rc.finish();
}

void BandFile::read_chunk(size_t c, vector<vector<int64_t>>& rows) const {
    const size_t first = c * h.rows_per_chunk, last = std::min<size_t>(h.n_snps, first + h.rows_per_chunk);
    const uint64_t a = chunk_off[c], b = chunk_off[c + 1];
    vector<uint8_t> buf(b - a);
    if (!buf.empty()) pread_all(fd, buf.data(), buf.size(), a, path);
    bcode::Models models; models.reset();
    bcode::Decoder rc(buf.data(), buf.size());
    const bool exact = h.round_bits == 0;
    const uint64_t N = h.n_samples;
    rows.resize(last - first);
    for (size_t j = first; j < last; j++) {
        const size_t l = snps.lo[j];
        vector<int64_t>& row = rows[j - first];
        row.resize(j - l);
        for (size_t i = l; i < j; i++) {
            if (exact) row[i - l] = bcode::get(rc, models.at(scale_ctx(snps, i, j), j - i)) + expect_C(snps.sxr(i), snps.sxr(j), N);
            else row[i - l] = bcode::get(rc, models.at(scale_ctx(snps, i, j), j - i));
        }
    }
}
