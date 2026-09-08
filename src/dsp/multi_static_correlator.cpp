#include <hampr/dsp/multi_static_correlator.hpp>
#include <hampr/utils/fftw_wrapper.hpp>
#include <hampr/utils/math_utils.hpp>
#include <cmath>
#include <algorithm>

namespace hampr {

MultiStaticCorrelator::CrossCorrelationResult
MultiStaticCorrelator::cross_correlate_pair(
        const IQBuffer& ref_a, const IQBuffer& ref_b,
        double fs, int site_a, int site_b) {
    CrossCorrelationResult result;
    result.site_a = site_a;
    result.site_b = site_b;

    int N = std::min(static_cast<int>(ref_a.size()),
                      static_cast<int>(ref_b.size()));
    if (N < 2) {
        result.time_delay_samples = 0.0;
        result.correlation_peak = 0.0;
        return result;
    }

    array a(ref_a.begin(), ref_a.begin() + N);
    array b(ref_b.begin(), ref_b.begin() + N);

    fft(a);
    fft(b);

    IQBuffer cross_spectrum(N);
    for (int i = 0; i < N; ++i)
        cross_spectrum[i] = a[i] * std::conj(b[i]);

    ifft(cross_spectrum);

    result.cross_spectrum = cross_spectrum;

    int peak_idx = 0;
    double max_val = 0.0;
    for (int i = 0; i < N; ++i) {
        double val = std::abs(cross_spectrum[i]);
        if (val > max_val) {
            max_val = val;
            peak_idx = i;
        }
    }

    int lag;
    if (peak_idx > N / 2)
        lag = peak_idx - N;
    else
        lag = peak_idx;

    result.time_delay_samples = static_cast<double>(lag);
    result.correlation_peak = max_val / N;
    return result;
}

std::vector<MultiStaticCorrelator::CrossCorrelationResult>
MultiStaticCorrelator::cross_correlate(const SynchronizedBatch& batch) {
    std::vector<CrossCorrelationResult> results;
    if (batch.aligned_data.size() < 2) return results;

    for (size_t i = 0; i < batch.aligned_data.size(); ++i) {
        for (size_t j = i + 1; j < batch.aligned_data.size(); ++j) {
            if (batch.aligned_data[i].empty() || batch.aligned_data[j].empty())
                continue;
            results.push_back(cross_correlate_pair(
                batch.aligned_data[i][0], batch.aligned_data[j][0],
                batch.fs, static_cast<int>(i), static_cast<int>(j)));
        }
    }
    return results;
}

std::vector<double>
MultiStaticCorrelator::compute_bistatic_ranges(
        const std::vector<CrossCorrelationResult>& results, double fs) {
    std::vector<double> ranges;
    for (const auto& r : results) {
        double delay_sec = r.time_delay_samples / fs;
        double range_m = 299792458.0 * std::abs(delay_sec);
        ranges.push_back(range_m);
    }
    return ranges;
}

double MultiStaticCorrelator::cooperative_confidence(
        const std::vector<CrossCorrelationResult>& results,
        const std::vector<IQBuffer>& site_surveillance) {
    if (results.empty() || site_surveillance.empty())
        return 0.0;

    double sum_peak = 0.0;
    for (const auto& r : results)
        sum_peak += r.correlation_peak;
    double avg_peak = sum_peak / results.size();

    double sum_energy = 0.0;
    for (const auto& ch : site_surveillance) {
        double energy = 0.0;
        for (const auto& s : ch)
            energy += std::abs(s) * std::abs(s);
        sum_energy += energy;
    }

    if (sum_energy <= 0.0)
        return 0.0;

    return avg_peak / (1.0 + std::sqrt(sum_energy));
}

} // namespace hampr
