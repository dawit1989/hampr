#ifndef HAMPR_DSP_MULTI_STATIC_FUSION_HPP
#define HAMPR_DSP_MULTI_STATIC_FUSION_HPP

#include <hampr/core/multi_site_types.hpp>
#include <hampr/core/types.hpp>
#include <hampr/dsp/localization.hpp>
#include <hampr/dsp/multi_static_correlator.hpp>
#include <vector>
#include <memory>

namespace hampr {

class MultiStaticFusion {
public:
    struct FusionResult {
        double latitude = 0.0;
        double longitude = 0.0;
        double altitude = 0.0;
        double residual_error = 0.0;
        FusionQuality quality;
        std::vector<FDoAMeasurement> fdoa_measurements;
    };

    MultiStaticFusion();

    // Estimate FDoA between two reference signals
    double estimate_fdoa(const IQBuffer& ref_a, const IQBuffer& ref_b, double fs);

    // Hybrid TDoA/FDoA localization
    FusionResult hybrid_localize(const SynchronizedBatch& batch,
                                 const TargetTrack& track,
                                 double carrier_freq_hz = 1e9);

    // GDOP for a given geometry (sites in ENU, target in ENU)
    double compute_gdop(const std::vector<std::pair<double, double>>& site_enu,
                        const std::pair<double, double>& target_enu);

    // Iterative refinement of position using TDoA and FDoA measurements
    std::pair<double, double> iterative_refine(
        const std::pair<double, double>& initial,
        const std::vector<double>& tdoa_samples,
        const std::vector<double>& fdoa_hz,
        const std::vector<std::pair<double, double>>& site_enu,
        double fs,
        double carrier_freq_hz,
        int max_iter = 3);

private:
    std::unique_ptr<TDoALocalizer> localizer_;
    std::unique_ptr<MultiStaticCorrelator> correlator_;

    static constexpr double SPEED_OF_LIGHT = 299792458.0;
    static constexpr double EARTH_RADIUS_M = 6371000.0;
};

} // namespace hampr

#endif // HAMPR_DSP_MULTI_STATIC_FUSION_HPP
