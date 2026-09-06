#ifndef HAMPR_DSP_METRIC_EXTRACTOR_HPP
#define HAMPR_DSP_METRIC_EXTRACTOR_HPP

#include <hampr/core/types.hpp>
#include <hampr/utils/math_utils.hpp>
#include <vector>

namespace hampr {

class MetricExtractor {
public:
    MetricExtractor() = default;

    double extract_snr(RDMap& rd_matrix,
                       std::vector<int> target_rd,
                       std::vector<int>& win) const;
};

} // namespace hampr

#endif // HAMPR_DSP_METRIC_EXTRACTOR_HPP
