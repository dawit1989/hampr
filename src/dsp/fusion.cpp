#include <hampr/dsp/fusion.hpp>
#include <hampr/core/exception.hpp>
#include <cmath>
#include <algorithm>

namespace hampr {

static constexpr double DEG_TO_RAD = 0.017453292519943295;
static constexpr double RAD_TO_DEG = 57.29577951308232;
static constexpr double EARTH_RADIUS_M = 6371000.0;

static void wgs84_to_enu(double lat, double lon,
                         double ref_lat, double ref_lon,
                         double& east, double& north) {
    double dlat = (lat - ref_lat) * DEG_TO_RAD;
    double dlon = (lon - ref_lon) * DEG_TO_RAD;
    double ref_lat_rad = ref_lat * DEG_TO_RAD;

    east  = dlon * std::cos(ref_lat_rad) * EARTH_RADIUS_M;
    north = dlat * EARTH_RADIUS_M;
}

static void enu_to_wgs84(double east, double north,
                         double ref_lat, double ref_lon,
                         double& lat, double& lon) {
    double dlat = north / EARTH_RADIUS_M;
    double dlon = east / (EARTH_RADIUS_M * std::cos(ref_lat * DEG_TO_RAD));
    lat = ref_lat + dlat * RAD_TO_DEG;
    lon = ref_lon + dlon * RAD_TO_DEG;
}

double ResultFusion::combined_sinr(const std::vector<ProcessingResult>& results) {
    if (results.empty())
        return 0.0;

    double sum_linear = 0.0;
    for (const auto& r : results)
        sum_linear += std::pow(10.0, r.sinr / 10.0);
    double avg_linear = sum_linear / results.size();
    return 10.0 * std::log10(avg_linear);
}

ResultFusion::FusedTarget ResultFusion::fuse_bearings(
        const std::vector<ProcessingResult>& results,
        const std::vector<ReceiverInfo>& receivers) {

    FusedTarget target;

    if (results.empty() || receivers.empty())
        return target;

    std::vector<std::pair<double, double>> site_positions;
    std::vector<double> bearings_deg;

    for (size_t i = 0; i < results.size() && i < receivers.size(); ++i) {
        const auto& r = results[i];
        if (r.sinr <= 0.0)
            continue;

        site_positions.emplace_back(receivers[i].latitude, receivers[i].longitude);
        bearings_deg.push_back(r.azimuth_estimated);
        target.contributing_bearings.push_back(r.azimuth_estimated);
        target.contributing_sites.push_back(receivers[i].id);
    }

    if (site_positions.size() < 2)
        return target;

    auto [lat, lon] = intersect_bearings(site_positions, bearings_deg);
    target.latitude = lat;
    target.longitude = lon;

    double alt_sum = 0.0;
    int alt_count = 0;
    for (size_t i = 0; i < results.size() && i < receivers.size(); ++i) {
        alt_sum += receivers[i].altitude;
        alt_count++;
    }
    target.altitude = alt_sum / alt_count;

    if (bearings_deg.size() == 2) {
        double angle_diff = std::abs(bearings_deg[0] - bearings_deg[1]);
        if (angle_diff > 180.0) angle_diff = 360.0 - angle_diff;
        target.confidence = std::sin(angle_diff * DEG_TO_RAD);
    } else {
        target.confidence = 1.0;
    }

    return target;
}

std::pair<double, double> ResultFusion::intersect_bearings(
        const std::vector<std::pair<double, double>>& site_positions,
        const std::vector<double>& bearings_deg) {

    if (site_positions.size() < 2)
        return {0.0, 0.0};

    double ref_lat = site_positions[0].first;
    double ref_lon = site_positions[0].second;

    std::vector<std::pair<double, double>> enu_positions;
    for (const auto& [lat, lon] : site_positions) {
        double e, n;
        wgs84_to_enu(lat, lon, ref_lat, ref_lon, e, n);
        enu_positions.emplace_back(e, n);
    }

    int n = static_cast<int>(enu_positions.size());

    if (n >= 2) {
        if (n > 2) {
            double sum_sin2 = 0, sum_cos2 = 0, sum_sin_cos = 0;
            double sum_x = 0, sum_y = 0;
            for (int i = 0; i < n; ++i) {
                double br = bearings_deg[i] * DEG_TO_RAD;
                double sb = std::sin(br), cb = std::cos(br);
                double ei = enu_positions[i].first;
                double ni = enu_positions[i].second;
                sum_sin2 += sb * sb;
                sum_cos2 += cb * cb;
                sum_sin_cos += sb * cb;
                sum_x += sb * (ei * sb - ni * cb);
                sum_y += (-cb) * (ei * sb - ni * cb);
            }

            double det = sum_sin2 * sum_cos2 - sum_sin_cos * sum_sin_cos;
            double x, y;
            if (std::abs(det) > 1e-20) {
                x = (sum_cos2 * sum_x - sum_sin_cos * sum_y) / det;
                y = (sum_sin2 * sum_y - sum_sin_cos * sum_x) / det;
            } else {
                x = enu_positions[0].first;
                y = enu_positions[0].second;
            }

            double lat, lon;
            enu_to_wgs84(x, y, ref_lat, ref_lon, lat, lon);
            return {lat, lon};
        } else {
            double b0 = bearings_deg[0] * DEG_TO_RAD;
            double b1 = bearings_deg[1] * DEG_TO_RAD;

            double e0 = enu_positions[0].first;
            double n0 = enu_positions[0].second;
            double e1 = enu_positions[1].first;
            double n1 = enu_positions[1].second;

            double sin_b0 = std::sin(b0), cos_b0 = std::cos(b0);
            double sin_b1 = std::sin(b1), cos_b1 = std::cos(b1);

            double c0 = e0 * sin_b0 - n0 * cos_b0;
            double c1 = e1 * sin_b1 - n1 * cos_b1;

            double det = sin_b0 * (-cos_b1) - sin_b1 * (-cos_b0);
            double x, y;
            if (std::abs(det) > 1e-20) {
                x = (c0 * (-cos_b1) - c1 * (-cos_b0)) / det;
                y = (sin_b0 * c1 - sin_b1 * c0) / det;
            } else {
                x = (e0 + e1) / 2.0;
                y = (n0 + n1) / 2.0;
            }

            double lat, lon;
            enu_to_wgs84(x, y, ref_lat, ref_lon, lat, lon);
            return {lat, lon};
        }
    }

    return {ref_lat, ref_lon};
}

} // namespace hampr
