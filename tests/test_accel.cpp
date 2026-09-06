#include <hampr/accel/fft_backend.hpp>
#include <hampr/accel/linalg_backend.hpp>
#include <hampr/utils/fftw_wrapper.hpp>
#include <hampr/utils/math_utils.hpp>
#include <hampr/utils/linalg_utils.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <iostream>
#include <cmath>
#include <cassert>
#include <random>
#include <algorithm>

static double max_rel_error(const hampr::array& a, const hampr::array& b) {
    assert(a.size() == b.size());
    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double denom = std::abs(b[i]);
        double err = denom > 1e-15 ? std::abs(a[i] - b[i]) / denom : std::abs(a[i] - b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

static double max_rel_error_mat(const hampr::mat& a, const hampr::mat& b) {
    assert(a.size() == b.size());
    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double err = max_rel_error(a[i], b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

int main() {
    using namespace hampr;
    std::cout << "Running Phase 4 backend tests..." << std::endl;

    auto fft_backend = create_fft_backend();
    auto linalg_backend = create_linalg_backend();

    // Test 1: FFTBackend forward matches direct fft()
    {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        array original(64);
        for (int i = 0; i < 64; ++i)
            original[i] = complex(dist(rng), dist(rng));

        array backend_result = original;
        array direct_result = original;

        fft_backend->forward(backend_result);
        fft(direct_result);

        double err = max_rel_error(backend_result, direct_result);
        assert(err < 1e-12);
        std::cout << "  Test 1 PASSED: FFTBackend forward matches direct fft" << std::endl;
    }

    // Test 2: FFTBackend inverse matches direct ifft()
    {
        std::mt19937 rng(99);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        array original(128);
        for (int i = 0; i < 128; ++i)
            original[i] = complex(dist(rng), dist(rng));

        array backend_result = original;
        array direct_result = original;

        fft_backend->inverse(backend_result);
        ifft(direct_result);

        double err = max_rel_error(backend_result, direct_result);
        assert(err < 1e-12);
        std::cout << "  Test 2 PASSED: FFTBackend inverse matches direct ifft" << std::endl;
    }

    // Test 3: FFTBackend round-trip is identity
    {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        array original(32);
        for (int i = 0; i < 32; ++i)
            original[i] = complex(dist(rng), dist(rng));

        array round_trip = original;
        fft_backend->forward(round_trip);
        fft_backend->inverse(round_trip);

        double err = max_rel_error(round_trip, original);
        assert(err < 1e-10);
        std::cout << "  Test 3 PASSED: FFT round-trip is identity" << std::endl;
    }

    // Test 4: LinAlgBackend multiply matches direct matMult
    {
        std::mt19937 rng(123);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        mat A(4, array(4));
        mat B(4, array(4));
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                A[i][j] = complex(dist(rng), dist(rng));
                B[i][j] = complex(dist(rng), dist(rng));
            }
        }

        mat backend_result = linalg_backend->multiply(A, B);
        mat direct_result = matMult(A, B);

        double err = max_rel_error_mat(backend_result, direct_result);
        assert(err < 1e-12);
        std::cout << "  Test 4 PASSED: LinAlgBackend multiply matches matMult" << std::endl;
    }

    // Test 5: LinAlgBackend inverse matches direct inverse
    {
        std::mt19937 rng(456);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        mat A(3, array(3));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                A[i][j] = complex(dist(rng), dist(rng));

        mat backend_result = linalg_backend->inverse(A);
        mat direct_result = inverse(A);

        double err = max_rel_error_mat(backend_result, direct_result);
        assert(err < 1e-10);
        std::cout << "  Test 5 PASSED: LinAlgBackend inverse matches direct inverse" << std::endl;
    }

    // Test 6: LinAlgBackend transpose_conj matches direct transposeConj
    {
        std::mt19937 rng(789);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        mat A(3, array(5));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 5; ++j)
                A[i][j] = complex(dist(rng), dist(rng));

        mat backend_result = linalg_backend->transpose_conj(A);
        mat direct_result = transposeConj(A);

        double err = max_rel_error_mat(backend_result, direct_result);
        assert(err < 1e-12);
        std::cout << "  Test 6 PASSED: LinAlgBackend transpose_conj matches transposeConj" << std::endl;
    }

    // Test 7: LinAlgBackend eigen_decompose produces valid decomposition
    {
        mat A(2, array(2));
        A[0][0] = complex(3.0, 0.0);
        A[0][1] = complex(1.0, 0.0);
        A[1][0] = complex(1.0, 0.0);
        A[1][1] = complex(3.0, 0.0);

        auto [eigenvalues, eigenvectors] = linalg_backend->eigen_decompose(A);

        assert(eigenvalues.size() == 2);
        assert(eigenvectors.size() == 2);
        assert(eigenvectors[0].size() == 2);

        std::vector<double> sorted_vals;
        for (const auto& ev : eigenvalues)
            sorted_vals.push_back(ev.real());
        std::sort(sorted_vals.begin(), sorted_vals.end());
        assert(std::abs(sorted_vals[0] - 2.0) < 1e-10);
        assert(std::abs(sorted_vals[1] - 4.0) < 1e-10);
        std::cout << "  Test 7 PASSED: eigen_decompose produces correct eigenvalues (2.0, 4.0)" << std::endl;
    }

    // Test 8: FlatBuffer element access and dimensions
    {
        FlatBuffer fb(3, 4);
        assert(fb.rows() == 3);
        assert(fb.cols() == 4);
        assert(fb.size() == 12);
        fb(1, 2) = complex(5.0, 3.0);
        assert(std::abs(fb(1, 2) - complex(5.0, 3.0)) < 1e-15);
        std::cout << "  Test 8 PASSED: FlatBuffer element access" << std::endl;
    }

    // Test 9: FlatBuffer from_mat / to_mat round-trip
    {
        mat m(3, array(4));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j)
                m[i][j] = complex(i + j, i * j);

        FlatBuffer fb;
        fb.from_mat(m);
        mat recovered = fb.to_mat();

        double err = max_rel_error_mat(m, recovered);
        assert(err < 1e-15);
        std::cout << "  Test 9 PASSED: FlatBuffer from_mat/to_mat round-trip (err=" << err << ")" << std::endl;
    }

    // Test 10: FlatBuffer extract_row / extract_column
    {
        FlatBuffer fb(3, 4);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j)
                fb(i, j) = complex(i, j);

        array row = fb.extract_row(1);
        assert(row.size() == 4);
        assert(std::abs(row[2] - complex(1.0, 2.0)) < 1e-15);

        array col = fb.extract_column(2);
        assert(col.size() == 3);
        assert(std::abs(col[1] - complex(1.0, 2.0)) < 1e-15);
        std::cout << "  Test 10 PASSED: FlatBuffer extract_row/extract_column" << std::endl;
    }

    // Test 11: LinAlgBackend flat buffer multiply matches mat multiply
    {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        mat A(4, array(4));
        mat B(4, array(4));
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                A[i][j] = complex(dist(rng), dist(rng));
                B[i][j] = complex(dist(rng), dist(rng));
            }
        }

        FlatBuffer fb_A; fb_A.from_mat(A);
        FlatBuffer fb_B; fb_B.from_mat(B);
        mat mat_result = matMult(A, B);

        FlatBuffer fb_result = linalg_backend->multiply_fb(fb_A, fb_B);
        mat fb_as_mat = fb_result.to_mat();

        double err = max_rel_error_mat(mat_result, fb_as_mat);
        assert(err < 1e-12);
        std::cout << "  Test 11 PASSED: multiply_fb matches matMult (err=" << err << ")" << std::endl;
    }

    // Test 12: LinAlgBackend flat buffer inverse matches mat inverse
    {
        std::mt19937 rng(777);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        mat A(3, array(3));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                A[i][j] = complex(dist(rng), dist(rng));

        FlatBuffer fb_A; fb_A.from_mat(A);
        mat mat_result = inverse(A);

        FlatBuffer fb_result = linalg_backend->inverse_fb(fb_A);
        mat fb_as_mat = fb_result.to_mat();

        double err = max_rel_error_mat(mat_result, fb_as_mat);
        assert(err < 1e-10);
        std::cout << "  Test 12 PASSED: inverse_fb matches inverse (err=" << err << ")" << std::endl;
    }

    // Test 13: LinAlgBackend flat buffer transpose_conj matches mat transpose_conj
    {
        std::mt19937 rng(333);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        mat A(3, array(5));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 5; ++j)
                A[i][j] = complex(dist(rng), dist(rng));

        FlatBuffer fb_A; fb_A.from_mat(A);
        mat mat_result = transposeConj(A);

        FlatBuffer fb_result = linalg_backend->transpose_conj_fb(fb_A);
        mat fb_as_mat = fb_result.to_mat();

        double err = max_rel_error_mat(mat_result, fb_as_mat);
        assert(err < 1e-12);
        std::cout << "  Test 13 PASSED: transpose_conj_fb matches transposeConj (err=" << err << ")" << std::endl;
    }

    // Test 14: FFTBackend flat buffer forward matches batch forward
    {
        std::mt19937 rng(55);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        int rows = 8, cols = 16;
        FlatBuffer fb(rows, cols);
        std::vector<complex> batch_data(rows * cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) {
                complex val(dist(rng), dist(rng));
                fb(i, j) = val;
                batch_data[i * cols + j] = val;
            }

        fft_backend->forward_fb(fb);
        fft_backend->forward_batch(batch_data, rows, cols);

        double max_err = 0.0;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) {
                double denom = std::abs(batch_data[i * cols + j]);
                double err = denom > 1e-15 ?
                    std::abs(fb(i, j) - batch_data[i * cols + j]) / denom :
                    std::abs(fb(i, j) - batch_data[i * cols + j]);
                if (err > max_err) max_err = err;
            }
        assert(max_err < 1e-12);
        std::cout << "  Test 14 PASSED: forward_fb matches forward_batch (err=" << max_err << ")" << std::endl;
    }

    std::cout << "All Phase 4 backend tests passed!" << std::endl;
    return 0;
}
