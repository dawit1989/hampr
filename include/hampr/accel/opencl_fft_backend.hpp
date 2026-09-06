#ifndef HAMPR_ACCEL_OPENCL_FFT_BACKEND_HPP
#define HAMPR_ACCEL_OPENCL_FFT_BACKEND_HPP

#include <hampr/accel/fft_backend.hpp>

namespace hampr {

// OpenCL fallback FFT backend.
// Used when a SYCL compiler is unavailable but OpenCL runtime is present.
class OpenCLFFTBackend : public FFTBackend {
public:
    void forward(array& data) override;
    void inverse(array& data) override;

    void forward_fb(FlatBuffer& data) override;
    void inverse_fb(FlatBuffer& data) override;
};

std::unique_ptr<FFTBackend> create_opencl_fft_backend();

} // namespace hampr
#endif // HAMPR_ACCEL_OPENCL_FFT_BACKEND_HPP
