#ifndef HAMPR_DSP_FUSION_HPP
#define HAMPR_DSP_FUSION_HPP

#include <hampr/core/multi_site_types.hpp>
#include <hampr/core/types.hpp>
#include <vector>
#include <string>

namespace hampr {

class ResultFusion {
public:
    struct FusedTarget {
        double latitude = 0.0;
        double longitude = 0.0;
        double altitude = 0.0;
        double confidence = 0.0;
        std::vector<double> contributing_bearings;
        std::vector<std::string> contributing_sites;
    };

    FusedTarget fuse_bearings(const std::vector<ProcessingResult>& results,
                               const std::vector<ReceiverInfo>& receivers);

    double combined_sinr(const std::vector<ProcessingResult>& results);

private:
    std::pair<double, double> intersect_bearings(
        const std::vector<std::pair<double, double>>& site_positions,
        const std::vector<double>& bearings_deg);

    static constexpr double EARTH_RADIUS_M = 6371000.0;
};

} // namespace hampr

#endif // HAMPR_DSP_FUSION_HPP
