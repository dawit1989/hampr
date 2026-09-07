#ifndef HAMPR_DSP_LOCALIZATION_HPP
#define HAMPR_DSP_LOCALIZATION_HPP

#include <hampr/core/multi_site_types.hpp>
#include <vector>
#include <utility>

namespace hampr {

class TDoALocalizer {
public:
    struct GeoResult {
        double latitude = 0.0;
        double longitude = 0.0;
        double altitude = 0.0;
        double residual_error = 0.0;
    };

    GeoResult localize(const SynchronizedBatch& batch,
                       const TargetTrack& track);

private:
    double estimate_tdoa(const IQBuffer& ref_a, const IQBuffer& ref_b, double fs);

    std::pair<double, double> tdoa_to_position(
        const std::vector<std::pair<double, double>>& site_positions,
        const std::vector<double>& tdoa_samples,
        double fs,
        double ref_lat, double ref_lon);

    static void wgs84_to_ecef(double lat, double lon, double alt,
                              double& x, double& y, double& z);
    static void ecef_to_wgs84(double x, double y, double z,
                              double& lat, double& lon, double& alt);

    static constexpr double SPEED_OF_LIGHT = 299792458.0;
    static constexpr double EARTH_RADIUS_M = 6371000.0;
};

} // namespace hampr

#endif // HAMPR_DSP_LOCALIZATION_HPP
