// compare_dumps.cpp - Compare two stage-dump trees produced by
// tests/numerical_dump.cpp (HAMPR) and validation/bench_dump.cpp (benchmark).
//
// Usage: compare_dumps <hampr_root> <bench_root>
// Both roots must contain a "dump/<idx>/<stage>.bin" tree.
// Writes a markdown report to stdout.
#include <algorithm>
#include <complex>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

struct Blob {
    int32_t rows = 0;
    int32_t cols = 0;
    std::vector<std::complex<double>> data;
};

static bool read_blob(const std::string& path, Blob& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char magic[4];
    if (!f.read(magic, 4)) return false;
    if (std::memcmp(magic, "HDMP", 4) != 0) {
        std::cerr << "bad magic in " << path << "\n";
        return false;
    }
    f.read(reinterpret_cast<char*>(&out.rows), 4);
    f.read(reinterpret_cast<char*>(&out.cols), 4);
    const std::size_t n = static_cast<std::size_t>(out.rows) * out.cols;
    out.data.resize(n);
    f.read(reinterpret_cast<char*>(out.data.data()), static_cast<std::streamsize>(n * sizeof(std::complex<double>)));
    return f.good() || f.eof();
}

// IEEE754 double ULP distance (treats +0/-0 as equal).
static int64_t ulp_diff(double a, double b) {
    int64_t ia, ib;
    std::memcpy(&ia, &a, sizeof(ia));
    std::memcpy(&ib, &b, sizeof(ib));
    if (a < 0) ia = 0x8000000000000000ULL - ia;
    else       ia = 0x8000000000000000ULL + ia;
    if (b < 0) ib = 0x8000000000000000ULL - ib;
    else       ib = 0x8000000000000000ULL + ib;
    return ia > ib ? ia - ib : ib - ia;
}

struct Stats {
    std::size_t total = 0;
    std::size_t exact = 0;
    double max_abs = 0.0;      // max |a-b| (complex magnitude)
    double max_abs_real = 0.0;
    double max_abs_imag = 0.0;
    double sum_sq = 0.0;
    double max_rel = 0.0;
    int64_t max_ulp = 0;
};

static Stats compare(const Blob& a, const Blob& b) {
    Stats s;
    if (a.rows != b.rows || a.cols != b.cols) {
        std::cerr << "  SHAPE MISMATCH " << a.rows << "x" << a.cols << " vs "
                  << b.rows << "x" << b.cols << "\n";
        return s;
    }
    s.total = a.data.size();
    for (std::size_t i = 0; i < s.total; ++i) {
        const double ra = a.data[i].real(), rb = b.data[i].real();
        const double ia = a.data[i].imag(), ib = b.data[i].imag();
        const double dr = ra - rb;
        const double di = ia - ib;
        const double mag = std::sqrt(dr * dr + di * di);
        if (mag == 0.0) s.exact++;
        if (mag > s.max_abs) s.max_abs = mag;
        if (std::fabs(dr) > s.max_abs_real) s.max_abs_real = std::fabs(dr);
        if (std::fabs(di) > s.max_abs_imag) s.max_abs_imag = std::fabs(di);
        if (mag > 0.0) {
            const double denom = std::sqrt(rb * rb + ib * ib);
            double rel = mag / (denom + 1e-300);
            if (rel > s.max_rel) s.max_rel = rel;
        }
        int64_t u = ulp_diff(ra, rb);
        if (u > s.max_ulp) s.max_ulp = u;
        u = ulp_diff(ia, ib);
        if (u > s.max_ulp) s.max_ulp = u;
        s.sum_sq += mag * mag;
    }
    return s;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: compare_dumps <hampr_root> <bench_root>\n";
        return 2;
    }
    const std::string hroot(argv[1]);
    const std::string broot(argv[2]);

    // Collect dataset indices.
    std::map<std::string, std::vector<std::string>> files_by_idx;
    {
        std::filesystem::path hdump = std::filesystem::path(hroot) / "dump";
        if (std::filesystem::exists(hdump)) {
            for (auto& e : std::filesystem::directory_iterator(hdump)) {
                if (e.is_directory())
                    for (auto& f : std::filesystem::directory_iterator(e.path()))
                        if (f.is_regular_file() && f.path().extension() == ".bin")
                            files_by_idx[e.path().filename().string()].push_back(f.path().filename().string());
            }
        }
    }
    for (auto& v : files_by_idx) std::sort(v.second.begin(), v.second.end());

    std::cout << "# Stage-by-Stage Numerical Comparison: HAMPR vs Benchmark\n\n";
    std::cout << "Comparison metric: max |HAMPR - benchmark| over all complex entries "
                 "(magnitude), plus exact-match count and IEEE754 double ULP distance.\n\n";

    struct Agg { std::size_t stages = 0; std::size_t all_exact = 0; };
    Agg agg;
    agg.all_exact = 1;

    for (const auto& [idx, stages] : files_by_idx) {
        std::cout << "## Dataset " << idx << "\n\n";
        std::cout << "| Stage | Rows x Cols | Elements | Exact | Max Abs Err | "
                     "Max |Re| | Max |Im| | RMS Err | Max Rel | Max ULP |\n";
        std::cout << "|-------|-------------|----------|-------|-------------|"
                     "---------|-----------------|---------|---------|\n";
        bool idx_all_exact = true;
        for (const auto& stage : stages) {
            std::string hp = (std::filesystem::path(hroot) / "dump" / idx / stage).string();
            std::string bp = (std::filesystem::path(broot) / "dump" / idx / stage).string();
            Blob ha, hb;
            if (!read_blob(hp, ha) || !read_blob(bp, hb)) {
                std::cout << "| " << stage << " | - | - | - | read failed |\n";
                continue;
            }
            Stats s = compare(ha, hb);
            const double rms = s.total ? std::sqrt(s.sum_sq / s.total) : 0.0;
            const bool exact = (s.exact == s.total && s.max_abs == 0.0);
            if (!exact) idx_all_exact = false;
            std::cout << std::setprecision(4);
            std::cout << "| " << stage << " | " << ha.rows << "x" << ha.cols
                      << " | " << s.total << " | " << s.exact << "/" << s.total
                      << " | " << std::scientific << std::setprecision(3) << s.max_abs
                      << " | " << s.max_abs_real << " | " << s.max_abs_imag
                      << " | " << std::fixed << std::setprecision(3) << rms
                      << " | " << std::scientific << std::setprecision(3) << s.max_rel
                      << " | " << s.max_ulp << " |\n";
        }
        std::cout << "\n**Dataset " << idx << " verdict:** "
                  << (idx_all_exact ? "EXACT (bit-identical at every stage)" : "DIVERGENCE DETECTED")
                  << "\n\n";
        agg.stages += stages.size();
        agg.all_exact &= idx_all_exact;
    }

    std::cout << "## Overall Verdict\n\n";
    std::cout << "**HAMPR is " << (agg.all_exact ? "BIT-IDENTICAL" : "NOT bit-identical")
              << " to the benchmark across all compared stages and datasets.**\n";
    return agg.all_exact ? 0 : 1;
}