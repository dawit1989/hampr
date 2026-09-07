#include <hampr/dsp/localization.hpp>
#include <hampr/core/exception.hpp>
#include <hampr/utils/fftw_wrapper.hpp>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace hampr {

static constexpr double DEG_TO_RAD = 0.017453292519943295;
static constexpr double RAD_TO_DEG = 57.29577951308232;
static constexpr double WGS84_A = 6378137.0;
static constexpr double WGS84_B = 6356752.314245;
static constexpr double WGS84_E2 = 1.0 - (WGS84_B * WGS84_B) / (WGS84_A * WGS84_A);

void TDoALocalizer::wgs84_to_ecef(double lat, double lon, double alt,
                                   double& x, double& y, double& z) {
    double lat_rad = lat * DEG_TO_RAD;
    double lon_rad = lon * DEG_TO_RAD;
    double N = WGS84_A / std::sqrt(1.0 - WGS84_E2 * std::sin(lat_rad) * std::sin(lat_rad));
    x = (N + alt) * std::cos(lat_rad) * std::cos(lon_rad);
    y = (N + alt) * std::cos(lat_rad) * std::sin(lon_rad);
    z = (N * (1.0 - WGS84_E2) + alt) * std::sin(lat_rad);
}

void TDoALocalizer::ecef_to_wgs84(double x, double y, double z,
                                   double& lat, double& lon, double& alt) {
    lon = std::atan2(y, x);
    double p = std::sqrt(x * x + y * y);
    double theta = std::atan2(z, p * (1.0 - WGS84_E2));
    double sin_theta = std::sin(theta);
    double N = WGS84_A / std::sqrt(1.0 - WGS84_E2 * sin_theta * sin_theta);

    lat = std::atan2(z + WGS84_E2 * N * sin_theta, p);
    double N_alt = WGS84_A / std::sqrt(1.0 - WGS84_E2 * std::sin(lat) * std::sin(lat));
    alt = p / std::cos(lat) - N_alt;

    lat *= RAD_TO_DEG;
    lon *= RAD_TO_DEG;
}

double TDoALocalizer::estimate_tdoa(const IQBuffer& ref_a, const IQBuffer& ref_b, double fs) {
    int N = std::min(static_cast<int>(ref_a.size()), static_cast<int>(ref_b.size()));
    if (N < 2)
        return 0.0;

    // Cross-correlation via FFT
    array a(ref_a.begin(), ref_a.begin() + N);
    array b(ref_b.begin(), ref_b.begin() + N);

    fft(a);
    fft(b);

    // Multiply A * conj(B)
    for (int i = 0; i < N; ++i)
        a[i] = a[i] * std::conj(b[i]);

    ifft(a);

    // Find peak
    int peak_idx = 0;
    double max_val = 0.0;
    for (int i = 0; i < N; ++i) {
        double val = std::abs(a[i]);
        if (val > max_val) {
            max_val = val;
            peak_idx = i;
        }
    }

    // Convert to signed lag
    int lag;
    if (peak_idx > N / 2)
        lag = peak_idx - N;
    else
        lag = peak_idx;

    return static_cast<double>(lag);
}

std::pair<double, double> TDoALocalizer::tdoa_to_position(
        const std::vector<std::pair<double, double>>& site_positions,
        const std::vector<double>& tdoa_samples,
        double fs,
        double ref_lat, double ref_lon) {

    if (site_positions.size() < 2)
        return {ref_lat, ref_lon};

    // Convert to ECEF
    std::vector<std::tuple<double, double, double>> ecef_pos;
    for (const auto& [lat, lon] : site_positions) {
        double x, y, z;
        wgs84_to_ecef(lat, lon, 0.0, x, y, z);
        ecef_pos.emplace_back(x, y, z);
    }

    // Set up least-squares for TDoA
    // |P - S_i| - |P - S_0| = c * tdoa_i / fs
    // Linearize around initial guess (reference site position)
    double x0, y0, z0;
    std::tie(x0, y0, z0) = ecef_pos[0];

    // Distance from initial guess to each site
    std::vector<double> d0;
    for (const auto& [x, y, z] : ecef_pos) {
        double dx = x0 - x;
        double dy = y0 - y;
        double dz = z0 - z;
        d0.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
    }

    // Build linearized system A * delta = b
    // Each equation: (x0 - sx)/d0_i * delta_x + ... - (x0 - s0x)/d0_0 * delta_x - ... = c*tdoa_i/fs + d0_i - d0_0
    int n_eq = static_cast<int>(tdoa_samples.size());
    std::vector<std::vector<double>> A(n_eq, std::vector<double>(3, 0.0));
    std::vector<double> b(n_eq, 0.0);

    double c_over_fs = SPEED_OF_LIGHT / fs;
    for (int i = 1; i < n_eq; ++i) {
        double si_x = std::get<0>(ecef_pos[i]);
        double si_y = std::get<1>(ecef_pos[i]);
        double si_z = std::get<2>(ecef_pos[i]);

        double di = d0[i] > 0 ? d0[i] : 1.0;
        double d0_ref = d0[0] > 0 ? d0[0] : 1.0;

        A[i-1][0] = (x0 - si_x) / di - (x0 - std::get<0>(ecef_pos[0])) / d0_ref;
        A[i-1][1] = (y0 - si_y) / di - (y0 - std::get<1>(ecef_pos[0])) / d0_ref;
        A[i-1][2] = (z0 - si_z) / di - (z0 - std::get<2>(ecef_pos[0])) / d0_ref;

        b[i-1] = c_over_fs * tdoa_samples[i] + di - d0_ref;
    }

    n_eq--; // we used n_eq-1 equations

    // Solve 3x3 normal equations: AtA * delta = Atb
    double AtA[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    double Atb[3] = {0, 0, 0};

    for (int i = 0; i < n_eq; ++i)
        for (int j = 0; j < 3; ++j) {
            Atb[j] += A[i][j] * b[i];
            for (int k = 0; k < 3; ++k)
                AtA[j][k] += A[i][j] * A[i][k];
        }

    // Solve via Cramer's rule
    double det = AtA[0][0] * (AtA[1][1] * AtA[2][2] - AtA[1][2] * AtA[2][1])
               - AtA[0][1] * (AtA[1][0] * AtA[2][2] - AtA[1][2] * AtA[2][0])
               + AtA[0][2] * (AtA[1][0] * AtA[2][1] - AtA[1][1] * AtA[2][0]);

    double delta_x, delta_y, delta_z;
    if (std::abs(det) > 1e-20) {
        double inv_det = 1.0 / det;
        delta_x = ((AtA[1][1] * AtA[2][2] - AtA[1][2] * AtA[2][1]) * Atb[0]
                 + (AtA[0][2] * AtA[2][1] - AtA[0][1] * AtA[2][2]) * Atb[1]
                 + (AtA[0][1] * AtA[1][2] - AtA[0][2] * AtA[1][1]) * Atb[2]) * inv_det;
        delta_y = ((AtA[1][2] * AtA[2][0] - AtA[1][0] * AtA[2][2]) * Atb[0]
                 + (AtA[0][0] * AtA[2][2] - AtA[0][2] * AtA[2][0]) * Atb[1]
                 + (AtA[0][2] * AtA[1][0] - AtA[0][0] * AtA[1][2]) * Atb[2]) * inv_det;
        delta_z = ((AtA[1][0] * AtA[2][1] - AtA[1][1] * AtA[2][0]) * Atb[0]
                 + (AtA[0][1] * AtA[2][0] - AtA[0][0] * AtA[2][1]) * Atb[1]
                 + (AtA[0][0] * AtA[1][1] - AtA[0][1] * AtA[1][0]) * Atb[2]) * inv_det;
    } else {
        delta_x = delta_y = delta_z = 0.0;
    }

    // Apply correction
    x0 += delta_x;
    y0 += delta_y;
    z0 += delta_z;

    double lat, lon, alt;
    ecef_to_wgs84(x0, y0, z0, lat, lon, alt);
    return {lat, lon};
}

TDoALocalizer::GeoResult TDoALocalizer::localize(const SynchronizedBatch& batch,
                                                   const TargetTrack& track) {
    GeoResult result;

    if (batch.aligned_data.size() < 2 || track.empty())
        return result;

    double fs = batch.fs;
    const auto& ref_track = track.front();

    // Estimate TDOA for each site pair (relative to site 0)
    std::vector<double> tdoa_samples;
    std::vector<std::pair<double, double>> site_positions;
    tdoa_samples.push_back(0.0);

    // Extract reference channel from each site
    std::vector<IQBuffer> ref_channels;
    for (const auto& site_data : batch.aligned_data) {
        if (!site_data.empty())
            ref_channels.push_back(site_data[0]);
        else
            ref_channels.push_back(IQBuffer());
    }

    for (size_t i = 1; i < ref_channels.size(); ++i) {
        double tdoa = estimate_tdoa(ref_channels[0], ref_channels[i], fs);
        tdoa_samples.push_back(tdoa);
    }

    for (const auto& info : batch.receivers)
        site_positions.emplace_back(info.latitude, info.longitude);

    double ref_lat = site_positions[0].first;
    double ref_lon = site_positions[0].second;

    auto [lat, lon] = tdoa_to_position(site_positions, tdoa_samples, fs, ref_lat, ref_lon);
    result.latitude = lat;
    result.longitude = lon;
    result.altitude = 0.0;

    // Residual error: distance from solution to each site vs expected
    double residual = 0.0;
    for (size_t i = 0; i < site_positions.size(); ++i) {
        double d_expected = SPEED_OF_LIGHT * tdoa_samples[i] / fs;
        double d_actual = 0.0;
        // Approximate with ENU distance
        // Simple Euclidean in lat/lon space
        double dlat = (lat - site_positions[i].first) * DEG_TO_RAD;
        double dlon = (lon - site_positions[i].second) * DEG_TO_RAD;
        double d_lat_m = dlat * EARTH_RADIUS_M;
        double d_lon_m = dlon * EARTH_RADIUS_M * std::cos(lat * DEG_TO_RAD);
        d_actual = std::sqrt(d_lat_m * d_lat_m + d_lon_m * d_lon_m);
        residual += std::abs(d_actual - std::abs(d_expected));
    }
    result.residual_error = residual / site_positions.size();

    return result;
}

} // namespace hampr

