#ifndef HAMPR_ACCEL_SYCL_LINALG_BACKEND_HPP
#define HAMPR_ACCEL_SYCL_LINALG_BACKEND_HPP

#include <hampr/accel/linalg_backend.hpp>

namespace hampr {

// SYCL-based linear algebra backend.
class SYCLLinAlgBackend : public LinAlgBackend {
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

std::unique_ptr<LinAlgBackend> create_sycl_linalg_backend();

} // namespace hampr
#endif // HAMPR_ACCEL_SYCL_LINALG_BACKEND_HPP
