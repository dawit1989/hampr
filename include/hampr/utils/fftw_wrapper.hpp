#ifndef HAMPR_UTILS_FFTW_WRAPPER_HPP
#define HAMPR_UTILS_FFTW_WRAPPER_HPP

#include <hampr/core/types.hpp>
#include <fftw3.h>
#include <cassert>
#include <unordered_map>

namespace hampr {

// Type aliases (complex, array) now centralized in hampr/core/types.hpp

namespace detail {

struct PlanKey {
    int n;
    int sign;
    bool operator==(const PlanKey& other) const { return n == other.n && sign == other.sign; }
};

struct PlanKeyHash {
    std::size_t operator()(const PlanKey& k) const {
        return std::hash<int>{}(k.n) ^ (std::hash<int>{}(k.sign) << 1);
    }
};

inline fftw_plan get_plan(int n, int sign) {
    static std::unordered_map<PlanKey, fftw_plan, PlanKeyHash> plan_cache;
    PlanKey key{n, sign};
    auto it = plan_cache.find(key);
    if (it != plan_cache.end()) {
        return it->second;
    }
    fftw_complex dummy[2];
    fftw_plan p = fftw_plan_dft_1d(n, dummy, dummy, sign, FFTW_ESTIMATE);
    plan_cache[key] = p;
    return p;
}

} // namespace detail

inline void fft(array& x) {
    int N = static_cast<int>(x.size());
    if (N == 0) return;
    fftw_complex* data = reinterpret_cast<fftw_complex*>(x.data());
    fftw_plan p = detail::get_plan(N, FFTW_FORWARD);
    fftw_execute_dft(p, data, data);
}

inline void ifft(array& x) {
    int N = static_cast<int>(x.size());
    if (N == 0) return;
    fftw_complex* data = reinterpret_cast<fftw_complex*>(x.data());
    fftw_plan p = detail::get_plan(N, FFTW_BACKWARD);
    fftw_execute_dft(p, data, data);
    complex scale(1.0 / N, 0.0);
    for (int i = 0; i < N; ++i)
        x[i] *= scale;
}

inline mat fftshift(const std::vector<array>& x) {
    int rows = static_cast<int>(x.size());
    int cols = static_cast<int>(x[0].size());
    std::vector<array> output = x;
    int row_shift = (rows + 1) / 2;
    int col_shift = (cols + 1) / 2;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            output[i][j] = x[(i + row_shift) % rows][(j + col_shift) % cols];
    return output;
}

} // namespace hampr

#endif // HAMPR_UTILS_FFTW_WRAPPER_HPP
