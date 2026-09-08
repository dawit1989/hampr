#ifndef HAMPR_DSP_MULTI_STATIC_CORRELATOR_HPP
#define HAMPR_DSP_MULTI_STATIC_CORRELATOR_HPP

#include <hampr/core/multi_site_types.hpp>
#include <hampr/core/types.hpp>
#include <vector>
#include <utility>

namespace hampr {

class MultiStaticCorrelator {
public:
    struct CrossCorrelationResult {
        int site_a;
        int site_b;
        double time_delay_samples;
        double correlation_peak;
        IQBuffer cross_spectrum;
    };

    std::vector<CrossCorrelationResult> cross_correlate(
        const SynchronizedBatch& batch);

    std::vector<double> compute_bistatic_ranges(
        const std::vector<CrossCorrelationResult>& results,
        double fs);

    double cooperative_confidence(
        const std::vector<CrossCorrelationResult>& results,
        const std::vector<IQBuffer>& site_surveillance);

private:
    CrossCorrelationResult cross_correlate_pair(
        const IQBuffer& ref_a,
        const IQBuffer& ref_b,
        double fs,
        int site_a, int site_b);
};

} // namespace hampr

#endif // HAMPR_DSP_MULTI_STATIC_CORRELATOR_HPP
