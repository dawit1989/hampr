#ifndef HAMPR_CORE_FLAT_BUFFER_HPP
#define HAMPR_CORE_FLAT_BUFFER_HPP

#include <hampr/core/types.hpp>
#include <vector>
#include <complex>
#include <cassert>

namespace hampr {

// FlatBuffer: contiguous 2D complex array for GPU-compatible memory layout.
// Unlike mat = vector<vector<complex>>, data is stored in a single
// contiguous buffer with row-major stride, enabling zero-copy GPU transfer
// and better CPU cache locality.
class FlatBuffer {
public:
    FlatBuffer() = default;
    FlatBuffer(int rows, int cols)
        : rows_(rows), cols_(cols), data_(static_cast<size_t>(rows) * cols) {}

    int rows() const { return rows_; }
    int cols() const { return cols_; }
    bool empty() const { return data_.empty(); }
    size_t size() const { return data_.size(); }

    // Row-major element access
    complex& operator()(int row, int col) {
        assert(row >= 0 && row < rows_ && col >= 0 && col < cols_);
        return data_[static_cast<size_t>(row) * cols_ + col];
    }
    const complex& operator()(int row, int col) const {
        assert(row >= 0 && row < rows_ && col >= 0 && col < cols_);
        return data_[static_cast<size_t>(row) * cols_ + col];
    }

    // Raw pointer access for GPU transfer
    complex* data() { return data_.data(); }
    const complex* data() const { return data_.data(); }
    std::vector<complex>& raw() { return data_; }
    const std::vector<complex>& raw() const { return data_; }

    // Extract a column as a flat vector
    array extract_column(int col) const {
        array result(rows_);
        for (int i = 0; i < rows_; ++i)
            result[i] = (*this)(i, col);
        return result;
    }

    // Extract a row as a flat vector
    array extract_row(int row) const {
        array result(cols_);
        for (int j = 0; j < cols_; ++j)
            result[j] = (*this)(row, j);
        return result;
    }

    // Set from a mat (vector of vectors) — for migration from old data model
    void from_mat(const mat& src) {
        rows_ = static_cast<int>(src.size());
        cols_ = static_cast<int>(src[0].size());
        data_.resize(static_cast<size_t>(rows_) * cols_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < cols_; ++j)
                (*this)(i, j) = src[i][j];
    }

    // Convert to mat (vector of vectors) — for compatibility with existing code
    mat to_mat() const {
        mat result(rows_);
        for (int i = 0; i < rows_; ++i) {
            result[i].resize(cols_);
            for (int j = 0; j < cols_; ++j)
                result[i][j] = (*this)(i, j);
        }
        return result;
    }

private:
    int rows_ = 0;
    int cols_ = 0;
    std::vector<complex> data_;
};

} // namespace hampr

#endif // HAMPR_CORE_FLAT_BUFFER_HPP
