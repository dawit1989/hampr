#include <hampr/dsp/multi_static_fusion.hpp>
#include <hampr/utils/fftw_wrapper.hpp>
#include <hampr/core/exception.hpp>
#include <hampr/utils/math_utils.hpp>
#include <limits>
#include <cmath>
#include <algorithm>

namespace hampr {

static const double DEG_TO_RAD = 0.017453292519943295;
static const double RAD_TO_DEG = 57.29577951308232;

MultiStaticFusion::MultiStaticFusion()
    : localizer_(std::make_unique<TDoALocalizer>()),
      correlator_(std::make_unique<MultiStaticCorrelator>()) {}

double MultiStaticFusion::estimate_fdoa(const IQBuffer& ref_a, const IQBuffer& ref_b, double fs) {
    return correlator_->estimate_fdoa(ref_a, ref_b, fs);
}

std::pair<double, double> MultiStaticFusion::iterative_refine(
        const std::pair<double, double>& initial,
        const std::vector<double>& tdoa_samples,
        const std::vector<double>& fdoa_hz,
        const std::vector<std::pair<double, double>>& site_enu,
        double fs, double carrier_freq_hz, int max_iter) {

    if (site_enu.empty() || tdoa_samples.empty())
        return initial;

    double x = initial.first;
    double y = initial.second;
    double c = SPEED_OF_LIGHT;
    double lam = c / carrier_freq_hz;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Build linearized system: H * delta = b
        std::vector<std::vector<double>> H;
        std::vector<double> b;

        // Reference site distance
        double dx0 = x - site_enu[0].first;
        double dy0 = y - site_enu[0].second;
        double d0 = std::sqrt(dx0 * dx0 + dy0 * dy0);
        if (d0 < 1e-6) d0 = 1e-6;

        // Unit vector from site 0 to target
        double sx0 = dx0 / d0;
        double sy0 = dy0 / d0;

        for (size_t i = 0; i < site_enu.size(); ++i) {
            double dx = x - site_enu[i].first;
            double dy = y - site_enu[i].second;
            double d = std::sqrt(dx * dx + dy * dy);
            if (d < 1e-6) d = 1e-6;

            // TDoa equation: d_i - d_0 = c * tdoa_i / fs
            // Linearized: (dx_i/d_i - dx_0/d_0) * delta_x + (dy_i/d_i - dy_0/d_0) * delta_y = c*tdoa_i/fs + d_0 - d_i
            H.push_back(std::vector<double>(2));
            H.back()[0] = dx / d - dx0 / d0;
            H.back()[1] = dy / d - dy0 / d0;
            b.push_back(c * tdoa_samples[i] / fs + d0 - d);

            // FDoa equation (if available)
            if (i < fdoa_hz.size() && std::abs(fdoa_hz[i]) > 1e-10) {
                // Unit vector from site i to target
                // Unit vector from site i to target (reversed from site-to-target)
                double sxi = dx / d;
                double syi = dy / d;
                // Geometry matrix row: (s_hat_i - s_hat_0) / lambda
                H.push_back(std::vector<double>(2));
                H.back()[0] = (sxi - sx0) / lam;
                H.back()[1] = (syi - sy0) / lam;
                b.push_back(fdoa_hz[i] - 0.0);  // expected FDoA (simplified: no velocity model)
            }
        }

        // Solve H * delta = b via normal equations (2x2)
        double AtA[2][2] = {{0, 0}, {0, 0}};
        double Atb[2] = {0, 0};

        for (size_t i = 0; i < H.size(); ++i) {
            AtA[0][0] += H[i][0] * H[i][0];
            AtA[0][1] += H[i][0] * H[i][1];
            AtA[1][0] += H[i][1] * H[i][0];
            AtA[1][1] += H[i][1] * H[i][1];
            Atb[0] += H[i][0] * b[i];
            Atb[1] += H[i][1] * b[i];
        }

        double det = AtA[0][0] * AtA[1][1] - AtA[0][1] * AtA[1][0];
        if (std::abs(det) < 1e-20)
            break;

        double delta_x = (AtA[1][1] * Atb[0] - AtA[0][1] * Atb[1]) / det;
        double delta_y = (AtA[0][0] * Atb[1] - AtA[1][0] * Atb[0]) / det;

        x += delta_x;
        y += delta_y;

        // Check convergence
        if (std::abs(delta_x) < 1e-3 && std::abs(delta_y) < 1e-3)
            break;
    }

    return {x, y};
}

double MultiStaticFusion::compute_gdop(
        const std::vector<std::pair<double, double>>& site_enu,
        const std::pair<double, double>& target_enu) {

    if (site_enu.size() < 2)
        return 0.0;

    // Build geometry matrix H (TDoA only for GDOP)
    double tx = target_enu.first;
    double ty = target_enu.second;

    std::vector<std::vector<double>> H;
    double dx0 = tx - site_enu[0].first;
    double dy0 = ty - site_enu[0].second;
    double d0 = std::sqrt(dx0 * dx0 + dy0 * dy0);
    if (d0 < 1e-6) d0 = 1e-6;
    double sx0 = dx0 / d0;
    double sy0 = dy0 / d0;

    for (size_t i = 1; i < site_enu.size(); ++i) {
        double dx = tx - site_enu[i].first;
        double dy = ty - site_enu[i].second;
        double d = std::sqrt(dx * dx + dy * dy);
        if (d < 1e-6) d = 1e-6;
        double sxi = dx / d;
        double syi = dy / d;

        H.push_back(std::vector<double>(2));
        H.back()[0] = sxi - sx0;
        H.back()[1] = syi - sy0;
    }

    // H^T * H (2x2)
    double AtA[2][2] = {{0, 0}, {0, 0}};
    for (size_t i = 0; i < H.size(); ++i) {
        AtA[0][0] += H[i][0] * H[i][0];
        AtA[0][1] += H[i][0] * H[i][1];
        AtA[1][0] += H[i][1] * H[i][0];
        AtA[1][1] += H[i][1] * H[i][1];
    }

    double det = AtA[0][0] * AtA[1][1] - AtA[0][1] * AtA[1][0];
    if (std::abs(det) < 1e-20)
        return std::numeric_limits<double>::infinity();

    // (H^T H)^-1
    double inv_00 = AtA[1][1] / det;
    double inv_11 = AtA[0][0] / det;

    // GDOP = sqrt(trace((H^T H)^-1))
    return std::sqrt(inv_00 + inv_11);
}

MultiStaticFusion::FusionResult MultiStaticFusion::hybrid_localize(
        const SynchronizedBatch& batch,
        const TargetTrack& track,
        double carrier_freq_hz) {

    FusionResult result;

    if (batch.aligned_data.size() < 2 || track.empty())
        return result;

    double fs = batch.fs;

    // Get TDoA measurements
    auto ms_results = correlator_->cross_correlate(batch);

    // Get FDoA measurements
    std::vector<FDoAMeasurement> fdoa_meas;
    std::vector<double> fdoa_hz;
    for (size_t i = 1; i < batch.aligned_data.size(); ++i) {
        if (batch.aligned_data[0].empty() || batch.aligned_data[i].empty())
            continue;
        double fdoa = correlator_->estimate_fdoa(
            batch.aligned_data[0][0], batch.aligned_data[i][0], fs);
        FDoAMeasurement m;
        m.site_a = 0;
        m.site_b = static_cast<int>(i);
        m.frequency_offset = fdoa;
        m.confidence = std::abs(fdoa) > 1.0 ? 1.0 : std::abs(fdoa) / 10.0;
        fdoa_meas.push_back(m);
        fdoa_hz.push_back(fdoa);
    }
    result.fdoa_measurements = fdoa_meas;

    // Get TDoA-only solution from TDoALocalizer
    TDoALocalizer::GeoResult tdoa_result = localizer_->localize(batch, track);

    // Convert site positions to ENU relative to reference site
    double ref_lat = batch.receivers[0].latitude;
    double ref_lon = batch.receivers[0].longitude;
    const double EARTH_R = EARTH_RADIUS_M;

    std::vector<std::pair<double, double>> site_enu;
    for (const auto& info : batch.receivers) {
        double dlat = (info.latitude - ref_lat) * DEG_TO_RAD;
        double dlon = (info.longitude - ref_lon) * DEG_TO_RAD;
        double east = dlon * std::cos(ref_lat * DEG_TO_RAD) * EARTH_R;
        double north = dlat * EARTH_R;
        site_enu.emplace_back(east, north);
    }

    // Convert TDoA solution to ENU
    double tdoa_dlat = (tdoa_result.latitude - ref_lat) * DEG_TO_RAD;
    double tdoa_dlon = (tdoa_result.longitude - ref_lon) * DEG_TO_RAD;
    std::pair<double, double> initial_enu(
        tdoa_dlon * std::cos(ref_lat * DEG_TO_RAD) * EARTH_R,
        tdoa_dlat * EARTH_R
    );

    // Build TDoA samples from correlator results
    std::vector<double> tdoa_samples;
    tdoa_samples.push_back(0.0);
    for (const auto& m : ms_results) {
        if (m.site_b > 0)
            tdoa_samples.push_back(m.time_delay_samples);
    }
    // Ensure correct ordering
    if (tdoa_samples.size() < batch.aligned_data.size()) {
        tdoa_samples.resize(batch.aligned_data.size(), 0.0);
    }

    // Get FDoA in correct ordering
    std::vector<double> fdoa_ordered(batch.aligned_data.size(), 0.0);
    for (const auto& m : fdoa_meas) {
        if (m.site_b < static_cast<int>(fdoa_ordered.size()))
            fdoa_ordered[m.site_b] = m.frequency_offset;
    }

    // Iterative refinement
    auto refined = iterative_refine(initial_enu, tdoa_samples, fdoa_ordered,
                                    site_enu, fs, carrier_freq_hz, 3);

    // Convert back to WGS84
    double lat_rad = ref_lat * DEG_TO_RAD;
    double dlat = refined.second / EARTH_R;
    double dlon = refined.first / (EARTH_R * std::cos(lat_rad));
    result.latitude = ref_lat + dlat * RAD_TO_DEG;
    result.longitude = ref_lon + dlon * RAD_TO_DEG;
    result.altitude = 0.0;

    // Compute GDOP
    result.quality.gdop = compute_gdop(site_enu, refined);
    result.quality.hdop = result.quality.gdop;
    result.quality.vdop = result.quality.gdop;
    result.quality.fdoa_confidence = fdoa_meas.empty() ? 0.0 : 1.0;
    result.quality.iterations = 3;

    result.residual_error = tdoa_result.residual_error;

    return result;
}

} // namespace hampr
