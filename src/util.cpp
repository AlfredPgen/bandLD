#include <fstream>

#include "common.hpp"

[[noreturn]] void die(const string& m) {
    fflush(stdout);
    fprintf(stderr, "Error : %s\n\n", m.c_str());
    exit(1);
}
void info(const char* k, const string& v) { printf("- %-22s: %s\n", k, v.c_str()); fflush(stdout); }
void warn(const string& m) { fprintf(stderr, "Warning : %s\n", m.c_str()); }

vector<string> split(const string& s, const string& delims, bool skip_empty) {
    vector<string> out; string cur;
    for (char c : s) {
        if (delims.find(c) != string::npos) { if (!skip_empty || !cur.empty()) out.push_back(cur); cur.clear(); }
        else cur += c;
    }
    if (!skip_empty || !cur.empty()) out.push_back(cur);
    return out;
}
string trim_cr(string s) { while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back(); return s; }

string norm_chrom(string c) {
    if (c.size() > 3 && (c.compare(0, 3, "chr") == 0 || c.compare(0, 3, "CHR") == 0)) c = c.substr(3);
    size_t z = 0;
    while (z + 1 < c.size() && c[z] == '0') z++;
    return c.substr(z);
}

vector<string> read_id_file_first_col(const string& fn, int skip_lines) {
    std::ifstream f(fn);
    if (!f) die("cannot open file '" + fn + "'");
    string line; vector<string> ids; int ln = 0;
    while (getline(f, line)) {
        if (++ln <= skip_lines) continue;
        auto v = split(trim_cr(line), " \t", true);
        if (!v.empty()) ids.push_back(v[0]);
    }
    return ids;
}

// sums over observed samples, missing index list, and the hard-call scale d:
// hard_d > 0 iff every observed k is in {0, d, 2d}
void finish_var(Var& v) {
    const size_t n = v.k.size();
    const int32_t* k = v.k.data();
    const uint8_t* m = v.miss.data();
    if (!v.sums_done) {
        int64_t sx = 0, sxx = 0; uint32_t nm = 0, c0 = 0, c2 = 0;
        const int32_t two = int32_t(2 * v.den);
        for (size_t i = 0; i < n; i++) {
            int64_t x = k[i]; sx += x; sxx += x * x; nm += m[i];
            c0 += (k[i] == 0) & (m[i] == 0); c2 += (k[i] == two) & (m[i] == 0);
        }
        v.sx = sx; v.sxx = sxx; v.n_miss = nm; v.c0 = c0; v.c2 = c2;
    }
    const uint32_t nm = v.n_miss;
    v.miss_idx.clear();
    if (nm) for (size_t i = 0; i < n; i++) if (m[i]) v.miss_idx.push_back((uint32_t)i);
    if (v.hard) { v.hard_d = v.den; return; }
    // general case: at most two distinct non-zero values a < b with b == 2a (or a single value)
    int32_t a = 0, b = 0;
    for (size_t i = 0; i < n; i++) {
        int32_t x = k[i];
        if (x == 0 || x == a || x == b) continue;
        if (!a) a = x; else if (!b) b = x; else { v.hard_d = 0; return; }
    }
    if (a && b) { if (a > b) std::swap(a, b); v.hard_d = (b == 2 * a) ? a : 0; }
    else v.hard_d = a ? a : v.den;   // single non-zero value (or none): k in {0, a}
}
