#ifndef HAMPR_ACCEL_FFT_BACKEND_HPP
#define HAMPR_ACCEL_FFT_BACKEND_HPP

#include <hampr/core/types.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <vector>
#include <memory>

namespace hampr {

class FFTBackend {
public:
    virtual ~FFTBackend() = default;

    // 1D array interface (CPU compatibility)
    virtual void forward(array& data) = 0;
    virtual void inverse(array& data) = 0;

    virtual void forward_batch(std::vector<complex>& data, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            array row(data.begin() + i * cols, data.begin() + (i + 1) * cols);
            forward(row);
            for (int j = 0; j < cols; ++j)
                data[i * cols + j] = row[j];
        }
    }
    virtual void inverse_batch(std::vector<complex>& data, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            array row(data.begin() + i * cols, data.begin() + (i + 1) * cols);
            inverse(row);
            for (int j = 0; j < cols; ++j)
                data[i * cols + j] = row[j];
        }
    }

    // FlatBuffer interface (GPU-compatible)
    virtual void forward_fb(FlatBuffer& data) {
        for (int i = 0; i < data.rows(); ++i) {
            array row(data.extract_row(i));
            forward(row);
            for (int j = 0; j < data.cols(); ++j)
                data(i, j) = row[j];
        }
    }
    virtual void inverse_fb(FlatBuffer& data) {
        for (int i = 0; i < data.rows(); ++i) {
            array row(data.extract_row(i));
            inverse(row);
            for (int j = 0; j < data.cols(); ++j)
                data(i, j) = row[j];
        }
    }
};

class FFTWBackend : public FFTBackend {
public:
    void forward(array& data) override;
    void inverse(array& data) override;
};

std::unique_ptr<FFTBackend> create_fft_backend();

} // namespace hampr
#endif // HAMPR_ACCEL_FFT_BACKEND_HPP
