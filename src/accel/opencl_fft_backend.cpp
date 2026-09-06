// OpenCL fallback FFT backend — conditionally compiled when HAMPR_ENABLE_GPU=OPENCL
//
// OpenCL C++ binding-based FFT implementation.
// Used when a SYCL compiler is unavailable but OpenCL runtime is present.
// Runs on any GPU vendor's platform (NVIDIA, AMD, Intel).
//
// Build: cmake -DHAMPR_ENABLE_GPU=OPENCL (requires OpenCL headers and libOpenCL)

#include <hampr/accel/opencl_fft_backend.hpp>

#define CL_HPP_CL_1_2_REQUIRED
#define CL_HPP_COMPLEX_TO_STR(x) #x
#define CL_HPP_TARGET_OPENCL_VERSION 120

#include <CL/cl2.hpp>
#include <stdexcept>

namespace hampr {

static cl_ulong find_device() {
    std::vector<cl::Device> devices;
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    for (const auto& platform : platforms) {
        std::vector<cl::Device> dev;
        platform.getDevices(CL_DEVICE_TYPE_GPU, &dev);
        if (!dev.empty())
            return 0;
        devices.insert(devices.end(), dev.begin(), dev.end());
    }
    throw std::runtime_error("No OpenCL GPU device found");
}

static void opencl_fft(array& data, bool inverse) {
    int N = static_cast<int>(data.size());
    if (N <= 1) return;

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.empty())
        throw std::runtime_error("No OpenCL platforms found");

    std::vector<cl::Device> devices;
    for (const auto& platform : platforms) {
        std::vector<cl::Device> dev;
        platform.getDevices(CL_DEVICE_TYPE_GPU, &dev);
        if (!dev.empty()) {
            devices = dev;
            break;
        }
    }
    if (devices.empty())
        throw std::runtime_error("No OpenCL GPU device found");

    cl::Context context(devices);
    cl::CommandQueue queue(context, devices[0]);

    // Naive DFT kernel
    std::string kernel_source = R"CLC(
        __kernel void dft(__global const float2* input,
                          __global float2* output,
                          const int N,
                          const int inverse) {
            int k = get_global_id(0);
            if (k >= N) return;
            float2 result = {0.0f, 0.0f};
            float sign = inverse ? 1.0f : -1.0f;
            for (int n = 0; n < N; ++n) {
                float angle = sign * 2.0f * 3.14159265358979f * k * n / N;
                float2 w = {cos(angle), sin(angle)};
                float2 x = input[n];
                result.x += w.x * x.x - w.y * x.y;
                result.y += w.x * x.y + w.y * x.x;
            }
            if (inverse) {
                result.x /= N;
                result.y /= N;
            }
            output[k] = result;
        }
    )CLC";

    cl::Program program(context, kernel_source);
    try {
        program.build({devices[0]});
    } catch (const cl::Error&) {
        throw std::runtime_error("OpenCL kernel build failed");
    }

    cl::Kernel kernel(program, "dft");

    // Create buffers
    cl::Buffer input_buf(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                         sizeof(cl_float2) * N,
                         reinterpret_cast<cl_float2*>(data.data()));
    std::vector<cl_float2> output_host(N);
    cl::Buffer output_buf(context, CL_MEM_WRITE_ONLY,
                          sizeof(cl_float2) * N);

    kernel.setArg(0, input_buf);
    kernel.setArg(1, output_buf);
    kernel.setArg(2, N);
    kernel.setArg(3, inverse ? 1 : 0);

    queue.enqueueNDRangeKernel(kernel, cl::NullRange,
                               cl::NDRange(N), cl::NullRange);
    queue.finish();

    queue.enqueueReadBuffer(output_buf, CL_TRUE, 0,
                            sizeof(cl_float2) * N, output_host.data());

    for (int i = 0; i < N; ++i)
        data[i] = complex(output_host[i].s[0], output_host[i].s[1]);
}

void OpenCLFFTBackend::forward(array& data) {
    opencl_fft(data, false);
}

void OpenCLFFTBackend::inverse(array& data) {
    opencl_fft(data, true);
}

void OpenCLFFTBackend::forward_fb(FlatBuffer& data) {
    for (int i = 0; i < data.rows(); ++i) {
        array row = data.extract_row(i);
        opencl_fft(row, false);
        for (int j = 0; j < data.cols(); ++j)
            data(i, j) = row[j];
    }
}

void OpenCLFFTBackend::inverse_fb(FlatBuffer& data) {
    for (int i = 0; i < data.rows(); ++i) {
        array row = data.extract_row(i);
        opencl_fft(row, true);
        for (int j = 0; j < data.cols(); ++j)
            data(i, j) = row[j];
    }
}

std::unique_ptr<FFTBackend> create_opencl_fft_backend() {
    return std::make_unique<OpenCLFFTBackend>();
}

} // namespace hampr
