#ifndef HAMPR_UTILS_LINALG_UTILS_HPP
#define HAMPR_UTILS_LINALG_UTILS_HPP

#include <hampr/core/types.hpp>
#include <Eigen/Eigenvalues>
#include <cassert>

namespace hampr {

// Type aliases (complex, mat) now centralized in hampr/core/types.hpp

inline mat inverse(const mat& x) {
    int n = static_cast<int>(x.size());
    assert(n > 1 && x.size() == x[0].size());
    Eigen::MatrixXcd temp(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            temp.row(i)[j] = x[i][j];
    temp = temp.inverse();
    mat output(n, std::vector<complex>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            output[i][j] = temp.row(i)[j];
    return output;
}

inline complex determinant(const mat& x) {
    int n = static_cast<int>(x.size());
    assert(n > 0 && x.size() == x[0].size());
    Eigen::MatrixXcd temp(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            temp.row(i)[j] = x[i][j];
    return complex(temp.determinant());
}

} // namespace hampr

#endif // HAMPR_UTILS_LINALG_UTILS_HPP
