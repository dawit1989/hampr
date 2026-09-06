#include <hampr/accel/fft_backend.hpp>
#include <hampr/utils/fftw_wrapper.hpp>

#ifdef HAMPR_GPU_BACKEND
#if HAMPR_GPU_BACKEND == SYCL
#include <hampr/accel/sycl_fft_backend.hpp>
#elif HAMPR_GPU_BACKEND == OPENCL
#include <hampr/accel/opencl_fft_backend.hpp>
#endif
#endif

namespace hampr {

void FFTWBackend::forward(array& data) {
    fft(data);
}

void FFTWBackend::inverse(array& data) {
    ifft(data);
}

std::unique_ptr<FFTBackend> create_fft_backend() {
#ifdef HAMPR_GPU_BACKEND
#if HAMPR_GPU_BACKEND == SYCL
    return create_sycl_fft_backend();
#elif HAMPR_GPU_BACKEND == OPENCL
    return create_opencl_fft_backend();
#endif
#endif
    return std::make_unique<FFTWBackend>();
}

} // namespace hampr
