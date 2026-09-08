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
    result.fdoa_hz = 0.0;

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

    // Estimate FDoA from phase difference between halves
    result.fdoa_hz = estimate_fdoa(ref_a, ref_b, fs);

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

double MultiStaticCorrelator::estimate_fdoa(const IQBuffer& ref_a, const IQBuffer& ref_b, double fs) {
    int N = std::min(static_cast<int>(ref_a.size()), static_cast<int>(ref_b.size()));
    if (N < 4)
        return 0.0;

    int half = N / 2;
    array a1(ref_a.begin(), ref_a.begin() + half);
    array b1(ref_b.begin(), ref_b.begin() + half);
    array a2(ref_a.begin() + half, ref_a.begin() + N);
    array b2(ref_b.begin() + half, ref_b.begin() + N);

    fft(a1); fft(b1);
    fft(a2); fft(b2);

    // Find peak frequency bin in combined cross-spectrum
    array combined(half);
    for (int i = 0; i < half; ++i)
        combined[i] = a1[i] * std::conj(b1[i]);
    double max_val = 0.0;
    int peak_bin = 0;
    for (int i = 0; i < half; ++i) {
        double val = std::abs(combined[i]);
        if (val > max_val) { max_val = val; peak_bin = i; }
    }

    // Phase difference at peak bin
    complex x1 = a1[peak_bin] * std::conj(b1[peak_bin]);
    complex x2 = a2[peak_bin] * std::conj(b2[peak_bin]);
    double delta_phase = std::arg(x2) - std::arg(x1);

    // Unwrap
    while (delta_phase > M_PI) delta_phase -= 2.0 * M_PI;
    while (delta_phase < -M_PI) delta_phase += 2.0 * M_PI;

    // FDoA = delta_phase / (2 * pi * T_half)
    double T_half = static_cast<double>(half) / fs;
    // Negate: cross-spectrum phase decreases when target frequency increases
    return -delta_phase / (2.0 * M_PI * T_half);
}

} // namespace hampr
