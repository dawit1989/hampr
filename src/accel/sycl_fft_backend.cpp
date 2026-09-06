// SYCL FFT backend — conditionally compiled when HAMPR_ENABLE_GPU=SYCL
//
// Single-source C++ FFT implementation using SYCL.
// Uses DPC++/hipSYCL/ComputeCpp to target CUDA, HIP, Level Zero, or Metal.
//
// Build: cmake -DHAMPR_ENABLE_GPU=SYCL (requires SYCL compiler)

#include <hampr/accel/sycl_fft_backend.hpp>
#include <sycl/sycl.hpp>
#include <stdexcept>

namespace hampr {

static void sycl_fft(array& data, bool inverse) {
    int N = static_cast<int>(data.size());
    if (N <= 1) return;

    sycl::queue queue(sycl::default_selector_v);

    // Host-to-device transfer
    using ComplexType = sycl::float2;
    sycl::float2* host_data = reinterpret_cast<sycl::float2*>(data.data());

    sycl::buffer<float, 1> buf(
        reinterpret_cast<float*>(host_data), sycl::range<1>(N * 2));

    queue.submit([&](sycl::handler& h) {
        auto data_acc = buf.get_access<sycl::access::mode::read_write>(h);

        h.parallel_for(sycl::range<1>(N), [=](sycl::id<1> idx) {
            int k = static_cast<int>(idx[0]);
            sycl::float2 result = {0.0f, 0.0f};
            for (int n = 0; n < N; ++n) {
                float angle = (inverse ? 2.0f : -2.0f) * 3.14159265358979f *
                              static_cast<float>(k) * static_cast<float>(n) /
                              static_cast<float>(N);
                sycl::float2 w = {cosf(angle), sinf(angle)};
                sycl::float2 x = data_acc[sycl::id<1>(n)];
                result.x += w.x * x.x - w.y * x.y;
                result.y += w.x * x.y + w.y * x.x;
            }
            if (inverse) {
                result.x /= static_cast<float>(N);
                result.y /= static_cast<float>(N);
            }
            data_acc[sycl::id<1>(k)] = result;
        });
    });

    queue.wait();

    // Copy back to host
    {
        auto host_acc = buf.get_access<sycl::access::mode::read>();
        for (int i = 0; i < N * 2; ++i)
            reinterpret_cast<float*>(host_data)[i] = host_acc[i];
    }
}

void SYCLFFTBackend::forward(array& data) {
    sycl_fft(data, false);
}

void SYCLFFTBackend::inverse(array& data) {
    sycl_fft(data, true);
}

void SYCLFFTBackend::forward_fb(FlatBuffer& data) {
    for (int i = 0; i < data.rows(); ++i) {
        array row = data.extract_row(i);
        sycl_fft(row, false);
        for (int j = 0; j < data.cols(); ++j)
            data(i, j) = row[j];
    }
}

void SYCLFFTBackend::inverse_fb(FlatBuffer& data) {
    for (int i = 0; i < data.rows(); ++i) {
        array row = data.extract_row(i);
        sycl_fft(row, true);
        for (int j = 0; j < data.cols(); ++j)
            data(i, j) = row[j];
    }
}

std::unique_ptr<FFTBackend> create_sycl_fft_backend() {
    return std::make_unique<SYCLFFTBackend>();
}

} // namespace hampr
