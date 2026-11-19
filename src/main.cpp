// bandLD: chromosome-wide banded LD in one compact file, with region and SNP-list look-up.
#include "build.hpp"
#include "kernels.hpp"

static void usage() {
    puts("bandLD " BANDLD_VERSION ": banded LD (all SNP pairs within a window) for a whole chromosome in one file\n\n"
         "bandLD build --bgen FILE [--bgi FILE] --out PREFIX [options]\n"
         "    --bgen FILE        BGEN v1.2/1.3 file (dosages or hard calls, biallelic, diploid)\n"
         "    --bgi FILE         bgenix index (default: FILE.bgi)\n"
         "    --chr NAME         chromosome (needed only if the index holds several)\n"
         "    --window-bp N      store pairs up to N base pairs apart (default 1000000)\n"
         "    --out PREFIX       writes PREFIX.bld, PREFIX.snps.tsv, PREFIX.skipped.tsv\n"
         "    --threads N        default 1\n"
         "    --round-bits B     store r rounded to 2^-B (default 15 = 16-bit precision, step 3e-5)\n"
         "    --exact            store exact integer cross products instead: r comes back as the same double\n"
         "                       fastLD computes; the file is larger for dosage data\n"
         "    --sample FILE      .sample file (needed with --incl if the BGEN has no sample IDs)\n"
         "    --incl FILE        sample IDs to use, one per line (default: all)\n"
         "    --extract FILE     variant IDs to use, one per line (default: all)\n"
         "    --min-maf X        skip variants with a lower minor allele frequency\n"
         "    --min-bp N, --max-bp N   restrict to a position range\n\n"
         "bandLD query --band PREFIX.bld --out FILE[.gz] (--min-bp N --max-bp N | --snps FILE) [options]\n"
         "    --chr NAME         optional check against the file's chromosome\n"
         "    --snps FILE        variant IDs, one per line; or a table with header 'rsid position allele1 allele2'\n"
         "                       (then alleles are checked and the sign follows the table's allele2)\n"
         "    --precision D      significant digits (default 6)\n"
         "    --threads N        default 1\n"
         "    writes the matrix as text (gzipped if FILE ends in .gz) and FILE.snps.tsv\n\n"
         "bandLD info --band PREFIX.bld\n"
         "bandLD --version");
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    const string cmd = argv[1];
    if (cmd == "--help" || cmd == "-h") { usage(); return 0; }
    if (cmd == "--version") {
        const char* k; ldk::pick_dense(&k);
        printf("bandLD %s\nkernel on this CPU: %s\n", BANDLD_VERSION, k);
        return 0;
    }
    auto need = [&](int& i) -> string {
        if (i + 1 >= argc) die(string("missing value for ") + argv[i]);
        return argv[++i];
    };
    auto num = [&](const string& flag, const string& s) {
        char* e; double v = strtod(s.c_str(), &e);
        if (s.empty() || *e) die("invalid value '" + s + "' for " + flag);
        return v;
    };
    printf("\n|--------------------------------------|\n| bandLD %-30s|\n|  banded LD, one file per chromosome  |\n|--------------------------------------|\n\n", BANDLD_VERSION);
    if (cmd == "build") {
        BuildOpts o;
        for (int i = 2; i < argc; i++) {
            const string a = argv[i];
            if (a == "--bgen") o.bgen = need(i);
            else if (a == "--bgi") o.bgi = need(i);
            else if (a == "--sample") o.sample = need(i);
            else if (a == "--incl") o.incl = need(i);
            else if (a == "--extract") o.extract = need(i);
            else if (a == "--chr") o.chrom = need(i);
            else if (a == "--out") o.out = need(i);
            else if (a == "--window-bp") o.window_bp = (uint64_t)num(a, need(i));
            else if (a == "--min-bp") o.min_bp = (uint32_t)num(a, need(i));
            else if (a == "--max-bp") o.max_bp = (uint32_t)num(a, need(i));
            else if (a == "--round-bits") o.round_bits = (int)num(a, need(i));
            else if (a == "--exact") o.round_bits = 0;
            else if (a == "--threads") o.threads = std::max(1, (int)num(a, need(i)));
            else if (a == "--min-maf") o.min_maf = num(a, need(i));
            else die("unknown argument '" + a + "' for build");
        }
        if (o.bgen.empty() || o.out.empty()) die("build needs --bgen and --out");
        if (o.bgi.empty()) o.bgi = o.bgen + ".bgi";
        if (o.round_bits < 0 || o.round_bits > 40) die("--round-bits must be between 1 and 40 (use --exact for exact storage)");
        if (o.window_bp < 1) die("--window-bp must be positive");
        run_build(o);
    } else if (cmd == "query") {
        QueryOpts o; bool have_min = false, have_max = false;
        for (int i = 2; i < argc; i++) {
            const string a = argv[i];
            if (a == "--band") o.band = need(i);
            else if (a == "--out") o.out = need(i);
            else if (a == "--snps") o.snps = need(i);
            else if (a == "--chr") o.chrom = need(i);
            else if (a == "--min-bp" || a == "--min_bp") { o.min_bp = (uint32_t)num(a, need(i)); have_min = true; }
            else if (a == "--max-bp" || a == "--max_bp") { o.max_bp = (uint32_t)num(a, need(i)); have_max = true; }
            else if (a == "--precision") o.precision = (int)num(a, need(i));
            else if (a == "--threads") o.threads = std::max(1, (int)num(a, need(i)));
            else die("unknown argument '" + a + "' for query");
        }
        if (o.band.empty() || o.out.empty()) die("query needs --band and --out");
        if (have_min != have_max) die("give both --min-bp and --max-bp");
        o.have_range = have_min;
        if (o.have_range && !o.snps.empty()) die("give either a position range or --snps, not both");
        if (o.precision < 1 || o.precision > 17) die("--precision must be between 1 and 17");
        run_query(o);
    } else if (cmd == "info") {
        string band;
        for (int i = 2; i < argc; i++) { const string a = argv[i]; if (a == "--band") band = need(i); else die("unknown argument '" + a + "' for info"); }
        if (band.empty()) die("info needs --band");
        run_info(band);
    } else die("unknown command '" + cmd + "' (build, query, info)");
    return 0;
}
