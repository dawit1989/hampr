#ifndef HAMPR_ACCEL_SYCL_FFT_BACKEND_HPP
#define HAMPR_ACCEL_SYCL_FFT_BACKEND_HPP

#include <hampr/accel/fft_backend.hpp>

namespace hampr {

// SYCL-based FFT backend.
// Single-source C++ — compiles to NVIDIA (CUDA), AMD (HIP), Intel (Level Zero),
// or Apple (Metal) via different SYCL compilers (DPC++, hipSYCL, ComputeCpp).
class SYCLFFTBackend : public FFTBackend {
public:
    void forward(array& data) override;
    void inverse(array& data) override;

    void forward_fb(FlatBuffer& data) override;
    void inverse_fb(FlatBuffer& data) override;
};

std::unique_ptr<FFTBackend> create_sycl_fft_backend();

} // namespace hampr
#endif // HAMPR_ACCEL_SYCL_FFT_BACKEND_HPP
