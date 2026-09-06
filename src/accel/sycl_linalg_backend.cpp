// SYCL LinAlg backend — conditionally compiled when HAMPR_ENABLE_GPU=SYCL
//
// Linear algebra backend using SYCL for GPU-accelerated matrix operations.
// Uses DPC++/hipSYCL/ComputeCpp to target CUDA, HIP, Level Zero, or Metal.
//
// For small matrices (e.g., 3x3 eigendecomposition in MUSIC), the CPU
// EigenBackend is used as fallback since GPU transfer overhead dominates.

#include <hampr/accel/sycl_linalg_backend.hpp>
#include <sycl/sycl.hpp>

namespace hampr {

static mat sycl_multiply(const mat& A, const mat& B) {
    int rows = static_cast<int>(A.size());
    int cols = static_cast<int>(B[0].size());
    int mid = static_cast<int>(B.size());

    sycl::queue queue(sycl::default_selector_v);

    // Flatten matrices for GPU transfer
    std::vector<double> a_data(rows * mid * 2);
    std::vector<double> b_data(mid * cols * 2);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < mid; ++j) {
            a_data[(i * mid + j) * 2] = A[i][j].real();
            a_data[(i * mid + j) * 2 + 1] = A[i][j].imag();
        }
    for (int i = 0; i < mid; ++i)
        for (int j = 0; j < cols; ++j) {
            b_data[(i * cols + j) * 2] = B[i][j].real();
            b_data[(i * cols + j) * 2 + 1] = B[i][j].imag();
        }

    sycl::buffer<double, 1> a_buf(a_data.data(), sycl::range<1>(rows * mid * 2));
    sycl::buffer<double, 1> b_buf(b_data.data(), sycl::range<1>(mid * cols * 2));
    std::vector<double> c_data(rows * cols * 2, 0.0);
    sycl::buffer<double, 1> c_buf(c_data.data(), sycl::range<1>(rows * cols * 2));

    queue.submit([&](sycl::handler& h) {
        auto a_acc = a_buf.get_access<sycl::access::mode::read>(h);
        auto b_acc = b_buf.get_access<sycl::access::mode::read>(h);
        auto c_acc = c_buf.get_access<sycl::access::mode::write>(h);

        h.parallel_for(sycl::range<2>(rows, cols), [=](sycl::id<2> idx) {
            int i = static_cast<int>(idx[0]);
            int j = static_cast<int>(idx[1]);
            double sum_r = 0.0, sum_i = 0.0;
            for (int k = 0; k < mid; ++k) {
                double a_r = a_acc[(i * mid + k) * 2];
                double a_i = a_acc[(i * mid + k) * 2 + 1];
                double b_r = b_acc[(k * cols + j) * 2];
                double b_i = b_acc[(k * cols + j) * 2 + 1];
                sum_r += a_r * b_r - a_i * b_i;
                sum_i += a_r * b_i + a_i * b_r;
            }
            c_acc[(i * cols + j) * 2] = sum_r;
            c_acc[(i * cols + j) * 2 + 1] = sum_i;
        });
    });

    queue.wait();

    mat result(rows, array(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[i][j] = complex(c_data[(i * cols + j) * 2],
                                   c_data[(i * cols + j) * 2 + 1]);
    return result;
}

mat SYCLLinAlgBackend::multiply(const mat& A, const mat& B) {
    return sycl_multiply(A, B);
}

mat SYCLLinAlgBackend::inverse(const mat& A) {
    // Small-matrix inverse: fallback to CPU EigenBackend
    EigenBackend cpu;
    return cpu.inverse(A);
}

std::pair<array, mat> SYCLLinAlgBackend::eigen_decompose(const mat& A) {
    // Eigendecomposition: fallback to CPU EigenBackend
    EigenBackend cpu;
    return cpu.eigen_decompose(A);
}

mat SYCLLinAlgBackend::transpose_conj(const mat& A) {
    int rows = static_cast<int>(A.size());
    int cols = static_cast<int>(A[0].size());
    mat result(cols, array(rows));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[j][i] = std::conj(A[i][j]);
    return result;
}

// --- FlatBuffer interface ---

FlatBuffer SYCLLinAlgBackend::multiply_fb(const FlatBuffer& A, const FlatBuffer& B) {
    mat A_mat = A.to_mat();
    mat B_mat = B.to_mat();
    mat result = sycl_multiply(A_mat, B_mat);
    FlatBuffer fb;
    fb.from_mat(result);
    return fb;
}

FlatBuffer SYCLLinAlgBackend::inverse_fb(const FlatBuffer& A) {
    EigenBackend cpu;
    mat A_mat = A.to_mat();
    mat result = cpu.inverse(A_mat);
    FlatBuffer fb;
    fb.from_mat(result);
    return fb;
}

std::pair<array, FlatBuffer> SYCLLinAlgBackend::eigen_decompose_fb(const FlatBuffer& A) {
    EigenBackend cpu;
    mat A_mat = A.to_mat();
    auto [eigenvalues, eigenvectors] = cpu.eigen_decompose(A_mat);
    FlatBuffer fb;
    fb.from_mat(eigenvectors);
    return {eigenvalues, fb};
}

FlatBuffer SYCLLinAlgBackend::transpose_conj_fb(const FlatBuffer& A) {
    mat A_mat = A.to_mat();
    mat result = transpose_conj(A_mat);
    FlatBuffer fb;
    fb.from_mat(result);
    return fb;
}

std::unique_ptr<LinAlgBackend> create_sycl_linalg_backend() {
    return std::make_unique<SYCLLinAlgBackend>();
}

} // namespace hampr
