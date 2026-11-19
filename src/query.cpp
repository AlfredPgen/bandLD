// bandLD query: LD matrix for a position range or a SNP list from a band file.
#include <libdeflate.h>
#include <omp.h>

#include <charconv>
#include <fstream>
#include <unordered_map>

#include "bandfile.hpp"
#include "build.hpp"

static void write_matrix(const string& fn, const vector<double>& R, size_t n, int precision, int nt) {
    FILE* f = fopen(fn.c_str(), "wb");
    if (!f) die("cannot write '" + fn + "'");
    const bool gz = fn.size() > 3 && fn.compare(fn.size() - 3, 3, ".gz") == 0;
    vector<libdeflate_compressor*> comp(gz ? size_t(nt) : size_t(0), nullptr);
    for (auto& c : comp) if (!(c = libdeflate_alloc_compressor(6))) die("out of memory (gzip compressor)");
    const size_t ROWS = 256, width = 8 + size_t(precision) + 8;
    vector<string> buf(ROWS), zbuf(comp.size());
    for (size_t a = 0; a < n; a += ROWS) {
        const size_t b = std::min(n, a + ROWS);
#pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
        for (size_t i = a; i < b; i++) {
            string& s = buf[i - a];
            s.resize(n * width + 2);
            char* p = s.data(); char* e = p + s.size();
            const double* row = &R[i * n];
            for (size_t j = 0; j < n; j++) {
                if (i == j) *p++ = '1';
                else if (std::isnan(row[j])) { memcpy(p, "nan", 3); p += 3; }
                else p = std::to_chars(p, e, row[j], std::chars_format::general, precision).ptr;
                *p++ = ' ';
            }
            *p++ = '\n';
            s.resize(p - s.data());
            if (gz) {
                const size_t th = size_t(omp_get_thread_num());
                string& z = zbuf[th];
                z.resize(libdeflate_gzip_compress_bound(comp[th], s.size()));
                const size_t zn = libdeflate_gzip_compress(comp[th], s.data(), s.size(), z.data(), z.size());
                if (!zn) die("gzip compression failed for '" + fn + "'");
                s.assign(z.data(), zn);
            }
        }
        for (size_t i = a; i < b; i++)
            if (fwrite(buf[i - a].data(), 1, buf[i - a].size(), f) != buf[i - a].size()) die("write error on '" + fn + "'");
    }
    for (auto c : comp) libdeflate_free_compressor(c);
    if (fclose(f) != 0) die("write error on '" + fn + "'");
}

void run_info(const string& band) {
    BandFile bf; bf.open(band);
    const BandHeader& h = bf.h;
    info("Band file", band);
    info("Chromosome", h.chrom);
    info("Samples", std::to_string(h.n_samples));
    info("Variants", std::to_string(h.n_snps));
    if (h.n_snps) info("Positions", std::to_string(bf.snps.pos.front()) + " to " + std::to_string(bf.snps.pos.back()));
    info("Window", std::to_string(h.window_bp) + " bp");
    info("Precision", h.round_bits ? "r rounded to 2^-" + std::to_string(h.round_bits) : string("exact (integer cross products)"));
    info("Pairs", std::to_string(h.n_pairs));
    char b[96]; snprintf(b, sizeof b, "%.3f GB, %.3f bytes/pair", double(bf.fsize) / 1e9, h.n_pairs ? double(bf.fsize) / double(h.n_pairs) : 0.0);
    info("File size", b);
}

void run_query(const QueryOpts& o) {
    const auto T0 = Clock::now();
    const int nt = std::max(1, o.threads);
    BandFile bf; bf.open(o.band);
    const BandHeader& h = bf.h; const BandSnps& S = bf.snps;
    const size_t M = S.size();
    info("Band file", o.band + " (chromosome " + h.chrom + ", " + std::to_string(M) + " variants, window " + std::to_string(h.window_bp) + " bp)");
    if (!o.chrom.empty() && norm_chrom(o.chrom) != norm_chrom(h.chrom))
        die("--chr " + o.chrom + " but '" + o.band + "' holds chromosome " + h.chrom);

    vector<uint32_t> sel; vector<int8_t> sgn;          // variant index and sign for each output row
    if (!o.snps.empty()) {
        std::ifstream f(o.snps);
        if (!f) die("cannot open SNP list '" + o.snps + "'");
        std::unordered_multimap<string, uint32_t> by;
        by.reserve(M * 2);
        for (uint32_t j = 0; j < M; j++) by.emplace(S.rsid[j], j);
        string line; int ir = 0, ip = -1, i1 = -1, i2 = -1; bool first = true; size_t ln = 0; bool bad = false;
        while (getline(f, line)) {
            ln++;
            auto t = split(trim_cr(line), " \t", true);
            if (t.empty()) continue;
            if (first) {
                first = false;
                auto col = [&](const char* n) { auto it = std::find(t.begin(), t.end(), n); return it == t.end() ? -1 : int(it - t.begin()); };
                if (col("rsid") >= 0) {                       // header line: rsid [chromosome] position allele1 allele2
                    ir = col("rsid"); ip = col("position"); i1 = col("allele1"); i2 = col("allele2");
                    if ((ip < 0) != (i1 < 0) || (i1 < 0) != (i2 < 0)) die("SNP list '" + o.snps + "': give all of position, allele1, allele2 or none");
                    continue;
                }
            }
            const int need = std::max(std::max(ir, ip), std::max(i1, i2));
            if ((int)t.size() <= need) die("SNP list '" + o.snps + "' line " + std::to_string(ln) + " is too short");
            const string& id = t[ir];
            vector<uint32_t> cand;
            int8_t sg = 1;
            if (ip >= 0) {                                    // match on position and the unordered allele pair
                const unsigned long p = strtoul(t[ip].c_str(), nullptr, 10);
                auto lo = std::lower_bound(S.pos.begin(), S.pos.end(), (uint32_t)p);
                for (auto it = lo; it != S.pos.end() && *it == p; ++it) {
                    const uint32_t j = uint32_t(it - S.pos.begin());
                    if ((S.a1[j] == t[i1] && S.a2[j] == t[i2]) || (S.a1[j] == t[i2] && S.a2[j] == t[i1])) cand.push_back(j);
                }
                if (cand.size() > 1) { vector<uint32_t> c2; for (uint32_t j : cand) if (S.rsid[j] == id) c2.push_back(j); cand = c2; }
                if (cand.size() == 1 && S.a2[cand[0]] != t[i2]) sg = -1;      // the file counts the list's allele1: flip
            } else {
                auto rg = by.equal_range(id);
                for (auto it = rg.first; it != rg.second; ++it) cand.push_back(it->second);
            }
            if (cand.size() != 1) {
                fprintf(stderr, "Error : SNP '%s' (line %zu of '%s') %s in the band file\n", id.c_str(), ln, o.snps.c_str(),
                        cand.empty() ? "was not found" : "matches several variants (give position and alleles)");
                bad = true; continue;
            }
            sel.push_back(cand[0]); sgn.push_back(sg);
        }
        if (bad) { fputc('\n', stderr); exit(1); }
        info("SNP list", o.snps + " (" + std::to_string(sel.size()) + " variants)");
    } else {
        if (!o.have_range) die("give --min-bp/--max-bp or --snps");
        const size_t a = std::lower_bound(S.pos.begin(), S.pos.end(), o.min_bp) - S.pos.begin();
        const size_t b = std::upper_bound(S.pos.begin(), S.pos.end(), o.max_bp) - S.pos.begin();
        for (size_t j = a; j < b; j++) { sel.push_back((uint32_t)j); sgn.push_back(1); }
        info("Region", std::to_string(o.min_bp) + " to " + std::to_string(o.max_bp) + " (" + std::to_string(sel.size()) + " variants)");
    }
    const size_t n = sel.size();
    if (!n) die("no variant selected");
    {
        char b[64]; snprintf(b, sizeof b, "%.2f GB", double(n) * double(n) * 8.0 / 1e9);
        info("Matrix", std::to_string(n) + " x " + std::to_string(n) + " (" + b + " in memory)");
    }

    // sorted unique variant indices and where each sits in the output
    vector<int32_t> where(M, -1);
    for (size_t x = 0; x < n; x++) {
        if (where[sel[x]] >= 0) die("SNP '" + S.rsid[sel[x]] + "' is listed twice");
        where[sel[x]] = (int32_t)x;
    }
    vector<uint32_t> u(sel); std::sort(u.begin(), u.end());
    vector<double> R(n * n, std::nan(""));
    for (size_t x = 0; x < n; x++) R[x * n + x] = 1.0;

    vector<size_t> chunks;
    for (uint32_t j : u) { const size_t c = j / h.rows_per_chunk; if (chunks.empty() || chunks.back() != c) chunks.push_back(c); }
    const bool exact = h.round_bits == 0;
    const double inv = std::ldexp(1.0, -int(h.round_bits));
    uint64_t stored = 0;
#pragma omp parallel num_threads(nt) reduction(+ : stored)
    {
        vector<vector<int64_t>> rows;
#pragma omp for schedule(dynamic, 1)
        for (size_t ci = 0; ci < chunks.size(); ci++) {
            const size_t c = chunks[ci], first = c * h.rows_per_chunk, last = std::min<size_t>(M, first + h.rows_per_chunk);
            bf.read_chunk(c, rows);
            auto jb = std::lower_bound(u.begin(), u.end(), (uint32_t)first), je = std::lower_bound(u.begin(), u.end(), (uint32_t)last);
            for (auto jt = jb; jt != je; ++jt) {
                const size_t j = *jt, l = S.lo[j];
                const vector<int64_t>& row = rows[j - first];
                const size_t xj = (size_t)where[j];
                for (auto it = std::lower_bound(u.begin(), jt, (uint32_t)l); it != jt; ++it) {
                    const size_t i = *it, xi = (size_t)where[i];
                    double r = exact ? r_exact(row[i - l], S, i, j, h.n_samples) : double(row[i - l]) * inv;
                    if (sgn[xi] * sgn[xj] < 0) r = -r;
                    if (r == 0) r = 0.0;                      // no "-0"
                    R[xi * n + xj] = r; R[xj * n + xi] = r;
                    stored++;
                }
            }
        }
    }
    const uint64_t want = uint64_t(n) * (n - 1) / 2;
    if (stored < want)
        warn(std::to_string(want - stored) + " of " + std::to_string(want) + " pairs are further apart than the window (" +
             std::to_string(h.window_bp) + " bp) and are written as nan");

    write_matrix(o.out, R, n, o.precision, nt);
    const string fsn = o.out + ".snps.tsv";
    FILE* ft = fopen(fsn.c_str(), "w");
    if (!ft) die("cannot write '" + fsn + "'");
    fprintf(ft, "#chrom\tpos\trsid\tallele1\tallele2\tmaf\n");
    for (size_t x = 0; x < n; x++) {
        const size_t j = sel[x]; const bool fl = sgn[x] < 0;
        fprintf(ft, "%s\t%u\t%s\t%s\t%s\t%.6g\n", h.chrom, S.pos[j], S.rsid[j].c_str(), (fl ? S.a2[j] : S.a1[j]).c_str(),
                (fl ? S.a1[j] : S.a2[j]).c_str(), S.maf(j, h.n_samples));
    }
    if (fclose(ft) != 0) die("write error on '" + fsn + "'");
    info("LD matrix", o.out);
    info("SNP table", fsn + " (LD sign refers to allele2)");
    char b[64]; snprintf(b, sizeof b, "%.2f s", secs(T0));
    info("Time", b);
}
