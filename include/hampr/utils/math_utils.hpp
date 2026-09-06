#ifndef HAMPR_UTILS_MATH_UTILS_HPP
#define HAMPR_UTILS_MATH_UTILS_HPP

#include <hampr/core/types.hpp>
#include <cassert>
#include <cmath>

namespace hampr {

// Type aliases (complex, array, mat) now centralized in hampr/core/types.hpp

inline array conjugate(const array& vec) {
    array output(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        output[i] = std::conj(vec[i]);
    return output;
}

inline mat conjugate(const mat& matrix) {
    mat output(matrix.size(), array(matrix[0].size()));
    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix[0].size(); ++j)
            output[i][j] = std::conj(matrix[i][j]);
    return output;
}

inline complex sum(const mat& matrix) {
    complex s(0.0, 0.0);
    for (const auto& row : matrix)
        for (const auto& val : row)
            s += val;
    return s;
}

inline complex dot(const array& a, const array& b) {
    assert(a.size() == b.size());
    complex s = 0;
    for (size_t i = 0; i < a.size(); ++i)
        s += a[i] * b[i];
    return s;
}

inline array dot(const mat& matrix, const array& arr) {
    assert(matrix[0].size() == arr.size());
    array output(matrix.size());
    for (size_t i = 0; i < matrix.size(); ++i)
        output[i] = dot(matrix[i], arr);
    return output;
}

inline mat reshape(const array& vec, int rows, int cols) {
    assert(rows * cols == static_cast<int>(vec.size()));
    mat output(rows, array(cols));
    int k = 0;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            output[i][j] = vec[k++];
    return output;
}

inline array shift(const array& x, int offset) {
    int N = static_cast<int>(x.size());
    if (!offset) return x;
    if (offset > N) return array(N, complex(0, 0));
    array output(N, complex(0, 0));
    for (int idx = N - 1; idx >= offset; --idx)
        output[idx] = x[idx - offset];
    return output;
}

inline mat transpose(const mat& x) {
    mat out(x[0].size(), array(x.size()));
    for (size_t i = 0; i < x.size(); ++i)
        for (size_t j = 0; j < x[0].size(); ++j)
            out[j][i] = x[i][j];
    return out;
}

inline mat transposeConj(const mat& x) {
    mat out(x[0].size(), array(x.size()));
    for (size_t i = 0; i < x.size(); ++i)
        for (size_t j = 0; j < x[0].size(); ++j)
            out[j][i] = std::conj(x[i][j]);
    return out;
}

inline mat matMult(const mat& lhs, const mat& rhs) {
    int rows = static_cast<int>(lhs.size());
    int cols = static_cast<int>(rhs[0].size());
    int mid = static_cast<int>(rhs.size());
    assert(static_cast<int>(lhs[0].size()) == mid);
    mat out(rows, array(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            for (int k = 0; k < mid; ++k)
                out[i][j] += lhs[i][k] * rhs[k][j];
    return out;
}

inline mat multiply(const mat& lhs, const mat& rhs) {
    int rows = static_cast<int>(lhs.size());
    int cols = static_cast<int>(lhs[0].size());
    assert(static_cast<int>(rhs.size()) == rows);
    assert(static_cast<int>(rhs[0].size()) == cols);
    mat out(rows, array(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            out[i][j] = lhs[i][j] * rhs[i][j];
    return out;
}

inline mat add(const mat& x, const mat& y) {
    int rows = static_cast<int>(x.size());
    int cols = static_cast<int>(x[0].size());
    mat out(rows, array(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            out[i][j] = x[i][j] + y[i][j];
    return out;
}

inline array convolve(const array& x, const array& y) {
    int N1 = static_cast<int>(x.size());
    int N2 = static_cast<int>(y.size());
    array out(N1 + N2 - 1);
    for (int i = 0; i < N1 + N2 - 1; ++i) {
        out[i] = 0;
        for (int j = 0; j < N2; ++j)
            if (i - j >= 0 && i - j < N1)
                out[i] += x[i - j] * y[j];
    }
    return out;
}

} // namespace hampr

#endif // HAMPR_UTILS_MATH_UTILS_HPP
