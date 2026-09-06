#ifndef HAMPR_ACCEL_LINALG_BACKEND_HPP
#define HAMPR_ACCEL_LINALG_BACKEND_HPP

#include <hampr/core/types.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <vector>
#include <memory>
#include <utility>

namespace hampr {

class LinAlgBackend {
public:
    virtual ~LinAlgBackend() = default;

    // mat/array interface (CPU compatibility)
    virtual mat multiply(const mat& A, const mat& B) = 0;
    virtual mat inverse(const mat& A) = 0;
    virtual std::pair<array, mat> eigen_decompose(const mat& A) = 0;
    virtual mat transpose_conj(const mat& A) = 0;

    // FlatBuffer interface (GPU-compatible)
    virtual FlatBuffer multiply_fb(const FlatBuffer& A, const FlatBuffer& B) {
        return FlatBuffer{};
    }
    virtual FlatBuffer inverse_fb(const FlatBuffer& A) {
        return FlatBuffer{};
    }
    virtual std::pair<array, FlatBuffer> eigen_decompose_fb(const FlatBuffer& A) {
        return {array{}, FlatBuffer{}};
    }
    virtual FlatBuffer transpose_conj_fb(const FlatBuffer& A) {
        return FlatBuffer{};
    }
};

class EigenBackend : public LinAlgBackend {
public:
    mat multiply(const mat& A, const mat& B) override;
    mat inverse(const mat& A) override;
    std::pair<array, mat> eigen_decompose(const mat& A) override;
    mat transpose_conj(const mat& A) override;

    FlatBuffer multiply_fb(const FlatBuffer& A, const FlatBuffer& B) override;
    FlatBuffer inverse_fb(const FlatBuffer& A) override;
    std::pair<array, FlatBuffer> eigen_decompose_fb(const FlatBuffer& A) override;
    FlatBuffer transpose_conj_fb(const FlatBuffer& A) override;
};

std::unique_ptr<LinAlgBackend> create_linalg_backend();

} // namespace hampr
#endif // HAMPR_ACCEL_LINALG_BACKEND_HPP
