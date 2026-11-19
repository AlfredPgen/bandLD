// bandLD build: banded exact LD for one chromosome, streamed along the chromosome.
//
// SNPs are decoded in position order in batches. Accepted SNPs are appended to a dense int16
// panel store (fastLD layout, kernels.hpp) that holds every SNP still within the window of the
// newest position. For each new SNP j the exact cross products with all stored earlier SNPs are
// computed by the integer kernels; the part inside the window becomes row j. Rows are entropy
// coded in chunks and written at once, so memory stays bounded by the window.
#include <omp.h>
#include <sys/mman.h>
#include <unistd.h>
#include <zstd.h>

#include <unordered_set>

#include "bandfile.hpp"
#include "build.hpp"
#include "kernels.hpp"

// common divisor of all numerators of a SNP (stops as soon as it reaches 1)
static uint32_t common_divisor(const vector<int32_t>& k) {
    uint32_t g = 0;
    for (int32_t x : k) {
        uint32_t a = (uint32_t)x;
        if (!a) continue;
        while (a) { const uint32_t t = g % a; g = a; a = t; }
        if (g == 1) break;
    }
    return g ? g : 1;
}

static string fmt_count(double v) {
    char b[64];
    if (v >= 1e12) snprintf(b, sizeof b, "%.2fT", v / 1e12);
    else if (v >= 1e9) snprintf(b, sizeof b, "%.2fB", v / 1e9);
    else if (v >= 1e6) snprintf(b, sizeof b, "%.2fM", v / 1e6);
    else if (v >= 1e3) snprintf(b, sizeof b, "%.1fK", v / 1e3);
    else snprintf(b, sizeof b, "%.0f", v);
    return b;
}
static string fmt_bytes(double v) {
    char b[64];
    if (v >= 1e12) snprintf(b, sizeof b, "%.2f TB", v / 1e12);
    else if (v >= 1e9) snprintf(b, sizeof b, "%.2f GB", v / 1e9);
    else if (v >= 1e6) snprintf(b, sizeof b, "%.1f MB", v / 1e6);
    else snprintf(b, sizeof b, "%.0f KB", v / 1e3);
    return b;
}
static string fmt_time(double s) {
    char b[64];
    if (s < 0 || !std::isfinite(s)) return "?";
    const long t = (long)(s + 0.5);
    if (t >= 86400) snprintf(b, sizeof b, "%ldd%02ldh", t / 86400, (t % 86400) / 3600);
    else if (t >= 3600) snprintf(b, sizeof b, "%ldh%02ldm", t / 3600, (t % 3600) / 60);
    else if (t >= 60) snprintf(b, sizeof b, "%ldm%02lds", t / 60, t % 60);
    else snprintf(b, sizeof b, "%lds", t);
    return b;
}

// dense panel store with a sliding front (see kernels.hpp for the layout)
struct Store {
    size_t n = 0, npair = 0, pstride = 0, cap_slots = 0, nslots = 0;
    int16_t* G = nullptr; size_t G_bytes = 0;
    void init(size_t n_, size_t cap) {
        n = n_; npair = (n + 1) / 2; pstride = 32 * std::max<size_t>(npair, 1);
        cap_slots = (cap + 15) / 16 * 16;
        G_bytes = std::max<size_t>(cap_slots / 16 * pstride * sizeof(int16_t), 4096);
        void* p = mmap(nullptr, G_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (p == MAP_FAILED) die("out of memory allocating the dosage store (" + fmt_bytes(double(G_bytes)) + ")");
        madvise(p, G_bytes, MADV_HUGEPAGE);
        G = (int16_t*)p;
    }
    ~Store() { if (G) munmap(G, G_bytes); }
    int16_t* panel(size_t p) const { return G + p * pstride; }
    size_t panel_bytes() const { return pstride * sizeof(int16_t); }
    // write the numerators of slot s divided by g (n values, each fits int16)
    void place(size_t s, const int32_t* k, int32_t g) {
        uint32_t* W = reinterpret_cast<uint32_t*>(panel(s / 16)) + (s % 16);
        const size_t np = n / 2;
        if (g == 1) for (size_t q = 0; q < np; q++) W[q * 16] = uint32_t(uint16_t(k[2 * q])) | (uint32_t(uint16_t(k[2 * q + 1])) << 16);
        else for (size_t q = 0; q < np; q++) W[q * 16] = uint32_t(uint16_t(k[2 * q] / g)) | (uint32_t(uint16_t(k[2 * q + 1] / g)) << 16);
        if (n & 1) W[np * 16] = uint32_t(uint16_t(k[n - 1] / g));
    }
    // drop the first d slots (d % 16 == 0): move the rest to the front and zero what is left behind
    void drop_front(size_t d) {
        if (!d) return;
        const size_t used = (nslots + 15) / 16, dp = d / 16;
        memmove(G, panel(dp), (used - dp) * panel_bytes());
        memset(panel(used - dp), 0, dp * panel_bytes());
        nslots -= d;
    }
};

void run_build(const BuildOpts& o) {
    const auto T0 = Clock::now();
    const int nt = std::max(1, o.threads);
    omp_set_num_threads(nt);

    // ---- variants
    string chrom; size_t n_multi = 0;
    vector<Snp> cand = bgi_variants(o.bgi, o.chrom, chrom, n_multi);
    info("BGEN file", o.bgen);
    info("BGI index", o.bgi);
    if (cand.empty()) die("no biallelic variant found in '" + o.bgi + "'" + (o.chrom.empty() ? "" : " on chromosome " + o.chrom));
    size_t n_listed = cand.size();
    if (!o.extract.empty()) {
        std::unordered_set<string> keep;
        for (auto& s : read_id_file_first_col(o.extract, 0)) keep.insert(s);
        vector<Snp> k2;
        for (auto& s : cand) if (keep.count(s.rsid)) k2.push_back(s);
        cand.swap(k2);
        info("Extract list", o.extract + " (" + std::to_string(keep.size()) + " IDs, " + std::to_string(cand.size()) + " found)");
        if (cand.empty()) die("no variant of the --extract list is in the index");
    }
    if (o.min_bp > 0 || o.max_bp > 0) {
        vector<Snp> k2;
        for (auto& s : cand) if (s.pos >= o.min_bp && (o.max_bp == 0 || s.pos <= o.max_bp)) k2.push_back(s);
        cand.swap(k2);
        if (cand.empty()) die("no variant in the requested position range");
    }
    const size_t Mc = cand.size();
    info("Chromosome", chrom);
    info("Variants in index", std::to_string(n_listed) + " biallelic" + (n_multi ? " (" + std::to_string(n_multi) + " multiallelic skipped)" : ""));
    info("Variants to process", std::to_string(Mc) + ", positions " + std::to_string(cand.front().pos) + " to " + std::to_string(cand.back().pos));

    std::unique_ptr<Source> src = open_bgen_list(o.bgen, o.sample, o.incl, cand);
    const size_t n = src->n_samples();
    info("Samples", std::to_string(n) + (n != src->n_file_samples() ? " of " + std::to_string(src->n_file_samples()) : ""));
    if (n < 2) die("at least 2 samples are needed");
    info("Window", std::to_string(o.window_bp) + " bp");
    info("Precision", o.round_bits ? "r rounded to 2^-" + std::to_string(o.round_bits) : string("exact (integer cross products)"));
    info("Threads", std::to_string(nt));

    // ---- size of the job from the candidate positions (upper bounds: filters come later)
    size_t maxwin = 0; double cand_pairs = 0;
    {
        size_t l = 0;
        for (size_t j = 0; j < Mc; j++) {
            while (uint64_t(cand[l].pos) + o.window_bp < uint64_t(cand[j].pos)) l++;
            maxwin = std::max(maxwin, j - l + 1); cand_pairs += double(j - l);
        }
    }
    const size_t BATCH = std::max<size_t>(32, std::min<size_t>(512, (size_t(256) << 20) / (5 * n)) / 16 * 16);
    Store st;
    st.init(n, maxwin + 2 * BATCH + 64);
    info("Largest window", std::to_string(maxwin) + " variants");
    info("Pairs (at most)", fmt_count(cand_pairs));
    info("Memory needed", "about " + fmt_bytes(double(st.G_bytes) + double(BATCH) * 5.0 * double(n) + double(BATCH) * 8.0 * double(maxwin) +
                                                 double(nt) * 256.0 * 8.0 * double(maxwin + 2 * BATCH) + 64e6));

    // ---- output
    const string fbld = o.out + ".bld", ftsv = o.out + ".snps.tsv", fskip = o.out + ".skipped.tsv";
    FILE* f = fopen(fbld.c_str(), "wb");
    if (!f) die("cannot write '" + fbld + "'");
    FILE* fs = fopen(fskip.c_str(), "w");
    if (!fs) die("cannot write '" + fskip + "'");
    fprintf(fs, "#chrom\tpos\trsid\tallele1\tallele2\treason\n");
    BandHeader h;
    h.round_bits = (uint32_t)o.round_bits; h.n_samples = n; h.window_bp = o.window_bp; h.rows_per_chunk = 64;
    snprintf(h.chrom, sizeof h.chrom, "%s", chrom.c_str());
    write_header(f, h);                                   // placeholder (table_off = 0 marks an unfinished file)
    info("Band file", fbld);
    printf("\n"); fflush(stdout);

    BandSnps S;                                           // accepted SNPs so far
    vector<uint64_t> chunk_off{HDR_BYTES};
    uint64_t file_pos = HDR_BYTES, n_pairs = 0;
    size_t g0 = 0;                                        // accepted index of store slot 0
    size_t rows_done = 0;                                 // rows already written
    vector<vector<int64_t>> pending;                      // rows rows_done .. (not yet written)
    size_t n_mono = 0, n_miss = 0, n_bits = 0, n_maf = 0, lpos = 0;
    double cand_done = 0;                                 // candidate pairs covered so far (for the ETA)

    const char* kname; ldk::DenseFn kern = ldk::pick_dense(&kname);
    vector<Var> vars(BATCH);
    vector<uint32_t> gdiv(BATCH, 1);
    vector<std::unique_ptr<ThreadCtx, void (*)(ThreadCtx*)>> ctx;
    for (int t = 0; t < nt; t++) ctx.push_back(src->make_ctx());
    vector<bcode::Models> models(nt);
    vector<vector<uint8_t>> cbuf;
    uint64_t maxk = 2;
    double t_dec = 0, t_ld = 0, t_enc = 0;
    auto t_last = Clock::now();
    const double every = getenv("BANDLD_PROGRESS_SEC") ? atof(getenv("BANDLD_PROGRESS_SEC")) : 30.0;
    const bool exact = o.round_bits == 0;
    const double rscale = std::ldexp(1.0, o.round_bits);

    auto flush_chunks = [&](bool final) {
        const size_t R = h.rows_per_chunk;
        const size_t nfull = final ? (pending.size() + R - 1) / R : pending.size() / R;
        if (!nfull) return;
        cbuf.resize(nfull);
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
        for (size_t c = 0; c < nfull; c++) {
            const size_t r0 = c * R, nr = std::min(R, pending.size() - r0);
            encode_chunk(h, S, rows_done + r0, pending, r0, nr, cbuf[c], models[omp_get_thread_num()]);
        }
        for (size_t c = 0; c < nfull; c++) {
            if (fwrite(cbuf[c].data(), 1, cbuf[c].size(), f) != cbuf[c].size()) die("write error on '" + fbld + "' (disk full?)");
            file_pos += cbuf[c].size(); chunk_off.push_back(file_pos);
        }
        const size_t nrows = std::min(pending.size(), nfull * R);
        pending.erase(pending.begin(), pending.begin() + nrows);
        rows_done += nrows;
    };

    for (size_t b0 = 0; b0 < Mc; b0 += BATCH) {
        const size_t b1 = std::min(Mc, b0 + BATCH);
        // -- decode
        auto tq = Clock::now();
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
        for (size_t c = b0; c < b1; c++) {
            src->decode(c, vars[c - b0], ctx[omp_get_thread_num()].get());
            gdiv[c - b0] = common_divisor(vars[c - b0].k);
        }
        t_dec += secs(tq); tq = Clock::now();

        // -- make room: slots that no later SNP can reach (all later positions are >= cand[b0].pos)
        {
            size_t d = 0;
            while (d < st.nslots && uint64_t(S.pos[g0 + d]) + o.window_bp < uint64_t(cand[b0].pos)) d++;
            d = d / 16 * 16;
            if (d >= 1024 || st.nslots + (b1 - b0) > st.cap_slots) { st.drop_front(d); g0 += d; }
            if (st.nslots + (b1 - b0) > st.cap_slots) die("internal error: window store too small");
        }
        // -- accept / reject, append to the store
        const size_t first_new = S.size();
        for (size_t c = b0; c < b1; c++) {
            const Var& v = vars[c - b0]; const Snp& sn = cand[c];
            const char* why = nullptr;
            const unsigned __int128 V = snp_V(v.sx, v.sxx, n);
            const double den = double(v.den), af = double(v.sx) / (2.0 * den * double(n)), maf = std::min(af, 1.0 - af);
            if (v.n_miss) { why = "missing genotypes"; n_miss++; }
            else if (v.bits > 10) { why = "probabilities with more than 10 bits"; n_bits++; }
            else if (V == 0) { why = "monomorphic"; n_mono++; }
            else if (maf < o.min_maf) { why = "MAF below --min-maf"; n_maf++; }
            if (why) { fprintf(fs, "%s\t%u\t%s\t%s\t%s\t%s\n", chrom.c_str(), sn.pos, sn.rsid.c_str(), sn.a1.c_str(), sn.a2.c_str(), why); continue; }
            const size_t j = S.size();
            const uint32_t g = gdiv[c - b0];
            S.pos.push_back(sn.pos); S.sx.push_back(v.sx); S.sxx.push_back(v.sxx); S.gcd.push_back(g); S.bits.push_back(v.bits);
            S.rsid.push_back(sn.rsid); S.a1.push_back(sn.a1); S.a2.push_back(sn.a2);
            S.vbits.push_back((uint8_t)bitlen128(V / ((unsigned __int128)g * g))); S.dV.push_back(double(V));
            while (uint64_t(S.pos[lpos]) + o.window_bp < uint64_t(sn.pos)) lpos++;
            S.lo.push_back(lpos);
            maxk = std::max<uint64_t>(maxk, 2 * uint64_t(v.den) / g);
            st.place(j - g0, v.k.data(), (int32_t)g);
            st.nslots = j - g0 + 1;
        }
        const size_t n_new = S.size() - first_new;
        // candidate pairs covered (for the ETA): pairs (i, c) of this batch within the window
        {
            static size_t cl = 0;
            for (size_t c = b0; c < b1; c++) { while (uint64_t(cand[cl].pos) + o.window_bp < uint64_t(cand[c].pos)) cl++; cand_done += double(c - cl); }
        }

        // -- cross products of the new slots with everything stored before them
        if (n_new) {
            const size_t s0 = first_new - g0, s1 = s0 + n_new;              // new slots [s0, s1)
            const size_t maxq = std::max<size_t>(1, (size_t(INT32_MAX) / (maxk * maxk)) / 2);
            const size_t QC = std::min<size_t>(256, maxq), IB = 48;
            // j-block width: at least two blocks per thread in a batch, 16..64 SNPs
            const size_t JB = std::max<size_t>(16, std::min<size_t>(64, (n_new / (2 * size_t(nt)) + 15) / 16 * 16));
            const size_t js0 = s0 / 16 * 16, nblk = (s1 - js0 + JB - 1) / JB;
            const size_t base = pending.size();
            pending.resize(base + n_new);
#pragma omp parallel num_threads(nt)
            {
                vector<int64_t> Lb;
#pragma omp for schedule(dynamic, 1)
                for (size_t blk = 0; blk < nblk; blk++) {
                    const size_t J0 = js0 + blk * JB, J1 = std::min(s1, J0 + JB), Ja = std::max(J0, s0);
                    if (Ja >= J1) continue;
                    const size_t J1r = (J1 + 7) / 8 * 8, nI = (J1 + 15) / 16, ldL = nI * 16;
                    const size_t ilo = S.lo[g0 + Ja] - g0;                   // smallest i-slot needed by this block
                    const size_t ib0 = ilo / 16;
                    Lb.assign((J1r - J0) * ldL, 0);
                    for (size_t q0 = 0; q0 < st.npair; q0 += QC) {
                        const size_t nq = std::min(QC, st.npair - q0);
                        for (size_t ib = ib0; ib < nI; ib += IB) {
                            const size_t ie = std::min(nI, ib + IB);
                            for (size_t js = J0; js < J1r; js += 8) {
                                const int16_t* B = st.panel(js / 16) + q0 * 32 + (js % 16) * 2;
                                for (size_t ig = ib; ig < ie; ig += 3) {
                                    if (ig * 16 >= js + 8) break;             // this and later groups are all above j
                                    const int IP = (int)std::min<size_t>(3, ie - ig);
                                    kern(IP, st.panel(ig) + q0 * 32, st.pstride, B, nq, &Lb[(js - J0) * ldL + ig * 16], ldL);
                                }
                            }
                        }
                    }
                    for (size_t j = Ja; j < J1; j++) {
                        const size_t gj = g0 + j, l = S.lo[gj] - g0;
                        vector<int64_t>& row = pending[base + (j - s0)];
                        row.resize(j - l);
                        const int64_t* L = &Lb[(j - J0) * ldL];
                        if (j == l) continue;                                 // first SNP of a window: empty row
                        if (exact) memcpy(row.data(), L + l, (j - l) * sizeof(int64_t));
                        else for (size_t i = l; i < j; i++) row[i - l] = (int64_t)std::llround(r_exact(L[i], S, g0 + i, gj, n) * rscale);
                    }
                }
            }
            for (size_t t = 0; t < n_new; t++) n_pairs += pending[base + t].size();
        }
        t_ld += secs(tq); tq = Clock::now();
        flush_chunks(false);
        t_enc += secs(tq);

        if (secs(t_last) >= every || b1 == Mc) {
            t_last = Clock::now();
            const double el = secs(T0), frac = cand_pairs > 0 ? cand_done / cand_pairs : double(b1) / double(Mc);
            const double eta = frac > 0 ? el * (1.0 - frac) / frac : -1;
            printf("[%5.1f%%] %zu / %zu variants | %s pairs | %s (%.2f bytes/pair) | %s pairs/s | elapsed %s | left %s\n",
                   100.0 * double(b1) / double(Mc), b1, Mc, fmt_count(double(n_pairs)).c_str(), fmt_bytes(double(file_pos)).c_str(),
                   n_pairs ? double(file_pos - HDR_BYTES) / double(n_pairs) : 0.0, fmt_count(el > 0 ? double(n_pairs) / el : 0).c_str(),
                   fmt_time(el).c_str(), b1 == Mc ? "0s" : fmt_time(eta).c_str());
            fflush(stdout);
        }
    }
    flush_chunks(true);
    fclose(fs);
    if (S.size() == 0) die("no variant was accepted (see '" + fskip + "')");

    // ---- SNP table, chunk index, header
    {
        const vector<uint8_t> ut = pack_snp_table(S);
        vector<uint8_t> ct(ZSTD_compressBound(ut.size()));
        const size_t cs = ZSTD_compress(ct.data(), ct.size(), ut.data(), ut.size(), 9);
        if (ZSTD_isError(cs)) die("zstd error on the SNP table");
        h.n_snps = S.size(); h.n_pairs = n_pairs; h.n_chunks = chunk_off.size() - 1;
        h.table_off = file_pos; h.table_csize = cs; h.table_usize = ut.size(); h.index_off = file_pos + cs;
        if (fwrite(ct.data(), 1, cs, f) != cs || fwrite(chunk_off.data(), 8, chunk_off.size(), f) != chunk_off.size())
            die("write error on '" + fbld + "' (disk full?)");
        file_pos += cs + chunk_off.size() * 8;
        write_header(f, h);
        if (fclose(f) != 0) die("write error on '" + fbld + "'");
    }
    // ---- SNP metadata table (for bgzip + tabix)
    {
        FILE* ft = fopen(ftsv.c_str(), "w");
        if (!ft) die("cannot write '" + ftsv + "'");
        fprintf(ft, "#chrom\tpos\trsid\tallele1\tallele2\tmaf\tindex\n");
        for (size_t j = 0; j < S.size(); j++)
            fprintf(ft, "%s\t%u\t%s\t%s\t%s\t%.6g\t%zu\n", chrom.c_str(), S.pos[j], S.rsid[j].c_str(), S.a1[j].c_str(), S.a2[j].c_str(), S.maf(j, n), j);
        if (fclose(ft) != 0) die("write error on '" + ftsv + "'");
    }
    printf("\n");
    info("Variants stored", std::to_string(S.size()));
    info("Variants skipped", std::to_string(n_mono + n_miss + n_bits + n_maf) + " (monomorphic " + std::to_string(n_mono) + ", missing data " +
                             std::to_string(n_miss) + ", >10-bit " + std::to_string(n_bits) + ", low MAF " + std::to_string(n_maf) + "): " + fskip);
    info("Pairs stored", fmt_count(double(n_pairs)) + " (" + std::to_string(n_pairs) + ")");
    info("Band file size", fmt_bytes(double(file_pos)) + (n_pairs ? " = " + std::to_string(double(file_pos) / double(n_pairs)).substr(0, 5) + " bytes/pair" : ""));
    info("Kernel", kname);
    {
        char b[160]; snprintf(b, sizeof b, "decode %s, LD %s, encode+write %s, total %s", fmt_time(t_dec).c_str(), fmt_time(t_ld).c_str(),
                              fmt_time(t_enc).c_str(), fmt_time(secs(T0)).c_str());
        info("Time", b);
    }
    info("SNP table", ftsv);
    printf("\nThe LD sign refers to allele2. To index the SNP table for region look-ups:\n  bgzip %s && tabix -s 1 -b 2 -e 2 %s.gz\n\n", ftsv.c_str(), ftsv.c_str());
}
