#ifndef HAMPR_DSP_SYNCHRONIZER_HPP
#define HAMPR_DSP_SYNCHRONIZER_HPP

#include <hampr/core/multi_site_types.hpp>
#include <hampr/accel/fft_backend.hpp>
#include <memory>

namespace hampr {

class Synchronizer {
public:
    explicit Synchronizer(double max_clock_error = 1e-4);

    SynchronizedBatch synchronize(const MultiSiteData& data);

    void calibrate(MultiSiteData& data);

    IQBuffer interpolate_shift(const IQBuffer& signal, double sample_offset);

private:
    double max_clock_error_;
    std::unique_ptr<FFTBackend> fft_backend_;
};

} // namespace hampr

#endif // HAMPR_DSP_SYNCHRONIZER_HPP
