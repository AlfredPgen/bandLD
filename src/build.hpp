#pragma once
#include "common.hpp"

struct BuildOpts {
    string bgen, bgi, sample, incl, chrom, out, extract;
    uint64_t window_bp = 1000000;
    uint32_t min_bp = 0, max_bp = 0;      // optional position range (0 = open)
    int round_bits = 15;                  // r rounded to 2^-15 (16-bit precision); 0 = exact
    int threads = 1;
    double min_maf = 0;
};
void run_build(const BuildOpts& o);

struct QueryOpts {
    string band, out, snps;               // snps: file with rsids, or a Z-style table (rsid chromosome position allele1 allele2)
    string chrom;
    bool have_range = false;
    uint32_t min_bp = 0, max_bp = 0;
    int precision = 6, threads = 1;
};
void run_query(const QueryOpts& o);
void run_info(const string& band);
