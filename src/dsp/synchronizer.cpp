#include <hampr/dsp/synchronizer.hpp>
#include <hampr/core/exception.hpp>
#include <cmath>
#include <limits>
#include <algorithm>

namespace hampr {

Synchronizer::Synchronizer(double max_clock_error)
    : max_clock_error_(max_clock_error),
      fft_backend_(create_fft_backend()) {}

void Synchronizer::calibrate(MultiSiteData& data) {
    for (size_t s = 0; s < data.iq_data.size(); ++s) {
        const ReceiverInfo& info = data.receivers[s];
        double amplitude = std::pow(10.0, info.gain / 20.0);
        complex phase_factor = complex(std::cos(info.phase_offset), std::sin(info.phase_offset));
        for (auto& ch : data.iq_data[s])
            for (auto& sample : ch)
                sample = sample * amplitude * phase_factor;
    }
}

IQBuffer Synchronizer::interpolate_shift(const IQBuffer& signal, double sample_offset) {
    int N = static_cast<int>(signal.size());
    if (std::abs(sample_offset) < 1e-12)
        return signal;

    int int_shift = static_cast<int>(std::floor(sample_offset));
    double frac = sample_offset - int_shift;

    IQBuffer shifted(N, complex(0, 0));

    for (int i = 0; i < N; ++i) {
        int src = i - int_shift;
        if (src < 0 || src >= N)
            continue;

        if (std::abs(frac) < 1e-6) {
            shifted[i] = signal[src];
        } else {
            int src_next = src + 1;
            if (src_next < N) {
                double w = 1.0 - frac;
                shifted[i] = signal[src] * w + signal[src_next] * frac;
            } else {
                shifted[i] = signal[src];
            }
        }
    }
    return shifted;
}

SynchronizedBatch Synchronizer::synchronize(const MultiSiteData& data) {
    if (data.iq_data.empty())
        return SynchronizedBatch{};

    SynchronizedBatch batch;
    batch.receivers = data.receivers;
    batch.fs = data.fs;
    batch.aligned_data.resize(data.iq_data.size());

    int min_samples = std::numeric_limits<int>::max();
    for (const auto& iq : data.iq_data)
        if (!iq.empty())
            min_samples = std::min(min_samples, static_cast<int>(iq[0].size()));
    batch.batch_samples = min_samples;

    for (size_t s = 0; s < data.iq_data.size(); ++s) {
        const ReceiverInfo& info = data.receivers[s];
        double sample_offset = info.clock_offset * data.fs;

        // Apply gain/phase calibration inline (avoids const_cast on input)
        double amplitude = std::pow(10.0, info.gain / 20.0);
        complex phase_factor = complex(std::cos(info.phase_offset), std::sin(info.phase_offset));
        complex cal_factor = amplitude * phase_factor;

        batch.aligned_data[s].resize(data.iq_data[s].size());
        for (size_t ch = 0; ch < data.iq_data[s].size(); ++ch) {
            const IQBuffer& channel = data.iq_data[s][ch];
            IQBuffer shifted = interpolate_shift(channel, sample_offset);

            // Apply calibration (linear op: order with interpolation is irrelevant)
            for (auto& sample : shifted)
                sample *= cal_factor;

            if (static_cast<int>(shifted.size()) > min_samples)
                shifted.resize(min_samples);
            else if (static_cast<int>(shifted.size()) < min_samples)
                shifted.resize(min_samples, complex(0, 0));

            batch.aligned_data[s][ch] = shifted;
        }
    }
    return batch;
}

} // namespace hampr
