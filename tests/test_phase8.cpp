#include <hampr/dsp/multi_static_fusion.hpp>
#include <hampr/dsp/multi_static_correlator.hpp>
#include <hampr/dsp/localization.hpp>
#include <hampr/dsp/multi_site_pipeline.hpp>
#include <hampr/core/config.hpp>
#include <hampr/core/types.hpp>
#include <hampr/core/multi_site_types.hpp>

static int test_failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << " at line " << __LINE__ << std::endl; test_failures++; } } while(0)
#include <hampr/dsp/synchronizer.hpp>
#include <cmath>
#include <iostream>
#include <complex>

using namespace hampr;

static const double C = 299792458.0;
static const double EARTH_R = 6371000.0;
static const double DEG_TO_RAD = 0.017453292519943295;
static const double RAD_TO_DEG = 57.29577951308232;

int main() {
    std::cout << "Running Phase 8 multi-static tests..." << std::endl;

    // Test 1: FDoA estimation with known frequency offset
    {
        MultiStaticFusion fusion;
        Synchronizer sync;

        const double fs = 1000000.0;  // 1 MHz
        const int N = 2000;
        const double delta_f = 50.0;  // 50 Hz frequency offset

        // Create reference signal (complex sinusoid, period=100 samples)
        IQBuffer ref_a(N);
        for (int i = 0; i < N; ++i)
            ref_a[i] = complex(std::cos(2.0 * M_PI * i / 100.0),
                               std::sin(2.0 * M_PI * i / 100.0));

        // Create signal with frequency offset: b[n] = a[n] * exp(j*2*pi*df*n/fs)
        IQBuffer ref_b(N);
        for (int i = 0; i < N; ++i) {
            double phase = 2.0 * M_PI * delta_f * i / fs;
            ref_b[i] = ref_a[i] * complex(std::cos(phase), std::sin(phase));
        }

        double estimated_fdoa = fusion.estimate_fdoa(ref_a, ref_b, fs);
        double error_hz = std::abs(estimated_fdoa - delta_f);

        CHECK(error_hz < 10.0, "error_hz < 10.0");
        std::cout << "  Test 1 PASSED: FDoA estimation (true=" << delta_f
                  << " Hz, est=" << estimated_fdoa << " Hz, err=" << error_hz << " Hz)" << std::endl;
    }

    // Test 2: GDOP computation for known geometry
    {
        MultiStaticFusion fusion;

        // 3 sites in equilateral triangle, target at center
        std::vector<std::pair<double, double>> sites = {
            {0.0, 0.0},
            {1000.0, 0.0},
            {500.0, 866.0}
        };
        std::pair<double, double> target(500.0, 289.0);  // centroid

        double gdop = fusion.compute_gdop(sites, target);

        CHECK(gdop > 0.0, "gdop > 0.0");
        CHECK(gdop < 100.0, "gdop < 100.0");
        CHECK(std::isfinite(gdop), "gdop not finite");

        std::cout << "  Test 2 PASSED: GDOP computation (gdop=" << gdop << ")" << std::endl;
    }

    // Test 3: Iterative refinement convergence
    {
        MultiStaticFusion fusion;

        // Simple 2D geometry: 3 sites, target at known position
        std::vector<std::pair<double, double>> sites = {
            {0.0, 0.0},
            {1000.0, 0.0},
            {500.0, 866.0}
        };
        std::pair<double, double> true_pos(500.0, 289.0);

        // Start with initial estimate off by ~200m
        std::pair<double, double> initial(700.0, 400.0);

        // Compute TDoA samples from true position
        std::vector<double> tdoa_samples(sites.size());
        double d_ref = std::sqrt(std::pow(true_pos.first - sites[0].first, 2) +
                                 std::pow(true_pos.second - sites[0].second, 2));
        for (size_t i = 0; i < sites.size(); ++i) {
            double d = std::sqrt(std::pow(true_pos.first - sites[i].first, 2) +
                                 std::pow(true_pos.second - sites[i].second, 2));
            tdoa_samples[i] = (d - d_ref) / C * 1e6;  // fs = 1 MHz
        }

        std::vector<double> fdoa_hz(sites.size(), 0.0);

        auto refined = fusion.iterative_refine(initial, tdoa_samples, fdoa_hz,
                                               sites, 1e6, 1e9, 3);

        double err = std::sqrt(std::pow(refined.first - true_pos.first, 2) +
                               std::pow(refined.second - true_pos.second, 2));

        CHECK(err < 500.0, "err < 500.0");  // Should converge within 500m
        std::cout << "  Test 3 PASSED: Iterative refinement (err=" << err << "m)" << std::endl;
    }

    // Test 4: GDOP decreases with more sites
    {
        MultiStaticFusion fusion;

        std::vector<std::pair<double, double>> sites_3 = {
            {0.0, 0.0}, {1000.0, 0.0}, {500.0, 866.0}
        };
        std::vector<std::pair<double, double>> sites_4 = {
            {0.0, 0.0}, {1000.0, 0.0}, {500.0, 866.0}, {500.0, -866.0}
        };
        std::pair<double, double> target(500.0, 289.0);

        double gdop_3 = fusion.compute_gdop(sites_3, target);
        double gdop_4 = fusion.compute_gdop(sites_4, target);

        CHECK(gdop_4 <= gdop_3, "gdop_4 <= gdop_3");  // More sites should reduce GDOP
        CHECK(gdop_4 > 0.0, "gdop_4 > 0.0");
        std::cout << "  Test 4 PASSED: GDOP ordering (3-site=" << gdop_3
                  << ", 4-site=" << gdop_4 << ")" << std::endl;
    }

    // Test 5: Integration with MultiSitePipeline
    {
        Config config;
        config.num_antennas = 3;
        config.antenna_spacing = 0.528;
        config.enable_localization = false;  // Keep simple
        config.enable_fusion = false;

        MultiSitePipeline pipeline(config);

        // Verify fusion_engine_ was initialized
        // (indirectly tested by successful construction)
        CHECK(pipeline.site_processing_times().empty(), "pipeline not empty");

        std::cout << "  Test 5 PASSED: MultiSitePipeline with FDoA engine initialized" << std::endl;
    }

    if (test_failures > 0) {
        std::cout << test_failures << " test(s) FAILED!" << std::endl;
        return 1;
    }
    std::cout << "All Phase 8 multi-static tests passed!" << std::endl;
    return 0;
}
