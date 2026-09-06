#include <hampr/accel/linalg_backend.hpp>
#include <hampr/utils/math_utils.hpp>
#include <hampr/utils/linalg_utils.hpp>
#include <Eigen/Eigenvalues>
#include <algorithm>

#ifdef HAMPR_GPU_BACKEND
#if HAMPR_GPU_BACKEND == SYCL
#include <hampr/accel/sycl_linalg_backend.hpp>
#endif
#endif

namespace hampr {

mat EigenBackend::multiply(const mat& A, const mat& B) {
    return matMult(A, B);
}

mat EigenBackend::inverse(const mat& A) {
    return hampr::inverse(A);
}

std::pair<array, mat> EigenBackend::eigen_decompose(const mat& A) {
    int n = static_cast<int>(A.size());
    Eigen::MatrixXcd mat_A(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat_A.row(i)[j] = A[i][j];

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(mat_A);

    array eigenvalues(n);
    mat eigenvectors(n, array(n));
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = ces.eigenvalues()[i];
        for (int j = 0; j < n; ++j)
            eigenvectors[j][i] = ces.eigenvectors().row(j)[i];
    }
    return {eigenvalues, eigenvectors};
}

mat EigenBackend::transpose_conj(const mat& A) {
    return hampr::transposeConj(A);
}

// --- FlatBuffer interface ---

FlatBuffer EigenBackend::multiply_fb(const FlatBuffer& A, const FlatBuffer& B) {
    mat A_mat = A.to_mat();
    mat B_mat = B.to_mat();
    mat result = matMult(A_mat, B_mat);
    FlatBuffer fb;
    fb.from_mat(result);
    return fb;
}

FlatBuffer EigenBackend::inverse_fb(const FlatBuffer& A) {
    mat A_mat = A.to_mat();
    mat result = hampr::inverse(A_mat);
    FlatBuffer fb;
    fb.from_mat(result);
    return fb;
}

std::pair<array, FlatBuffer> EigenBackend::eigen_decompose_fb(const FlatBuffer& A) {
    mat A_mat = A.to_mat();
    auto [eigenvalues, eigenvectors] = eigen_decompose(A_mat);
    FlatBuffer fb;
    fb.from_mat(eigenvectors);
    return {eigenvalues, fb};
}

FlatBuffer EigenBackend::transpose_conj_fb(const FlatBuffer& A) {
    mat A_mat = A.to_mat();
    mat result = hampr::transposeConj(A_mat);
    FlatBuffer fb;
    fb.from_mat(result);
    return fb;
}

std::unique_ptr<LinAlgBackend> create_linalg_backend() {
#ifdef HAMPR_GPU_BACKEND
#if HAMPR_GPU_BACKEND == SYCL
    return create_sycl_linalg_backend();
#endif
#endif
    return std::make_unique<EigenBackend>();
}

} // namespace hampr
