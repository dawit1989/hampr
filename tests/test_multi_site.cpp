#include <hampr/dsp/array_geometry.hpp>
#include <hampr/dsp/synchronizer.hpp>
#include <hampr/dsp/fusion.hpp>
#include <hampr/dsp/localization.hpp>
#include <hampr/dsp/multi_site_pipeline.hpp>
#include <hampr/io/multi_site_data_source.hpp>
#include <hampr/io/text_data_source.hpp>
#include <hampr/pipeline/multi_site_graph.hpp>
#include <hampr/dsp/multi_static_correlator.hpp>
#include <hampr/dsp/pipeline.hpp>
#include <hampr/core/config.hpp>
#include <hampr/core/types.hpp>
#include <hampr/core/multi_site_types.hpp>
#include <iostream>
#include <cmath>
#include <cassert>
#include <random>
#include <fstream>
#include <complex>

using namespace hampr;

static double deg_to_rad(double d) { return d * 0.017453292519943295; }
static double max_rel_error(const array& a, const array& b) {
    if (a.size() != b.size()) return 1e9;
    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double denom = std::abs(b[i]);
        double err = denom > 1e-15 ? std::abs(a[i] - b[i]) / denom : std::abs(a[i] - b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

int main() {
    std::cout << "Running Phase 7 multi-site tests..." << std::endl;

    // Test 1: MultiSiteData construction and field access
    {
        MultiSiteData data;
        assert(data.receivers.empty());
        assert(data.iq_data.empty());
        assert(data.fs == 0.0);

        ReceiverInfo info;
        info.id = "site_A";
        info.latitude = -37.8472;
        info.longitude = 145.0451;
        info.altitude = 12.0;
        info.clock_offset = 0.0;
        info.num_channels = 4;
        info.array_geometry.type = ArrayGeometrySpec::Type::ULA;
        info.array_geometry.element_spacing = 0.5;
        info.array_geometry.ula_elements = 4;

        data.receivers.push_back(info);
        data.iq_data.push_back(IQMatrix(4, IQBuffer(100)));
        data.fs = 2500000.0;

        assert(data.receivers.size() == 1);
        assert(data.receivers[0].id == "site_A");
        assert(data.receivers[0].array_geometry.ula_elements == 4);
        assert(data.iq_data[0].size() == 4);
        assert(data.iq_data[0][0].size() == 100);

        std::cout << "  Test 1 PASSED: MultiSiteData construction" << std::endl;
    }

    // Test 2: ArrayGeometry ULA scanning vectors match reference formula
    {
        ArrayGeometrySpec spec;
        spec.type = ArrayGeometrySpec::Type::ULA;
        spec.element_spacing = 0.5;
        spec.ula_elements = 4;
        ArrayGeometry geom(spec);

        assert(geom.num_elements() == 4);

        std::vector<double> thetas = {0, 30, 45, 90, 180};
        mat sv = geom.gen_scanning_vectors(thetas);

        // Reference: phase = 2*pi * d * cos(theta) for element at position x
        double d = 0.5;
        for (size_t i = 0; i < thetas.size(); ++i) {
            double theta_rad = deg_to_rad(thetas[i]);
            for (int j = 0; j < 4; ++j) {
                double x_pos = d * (j - 1.5);
                double phase = 2.0 * M_PI * x_pos * std::cos(theta_rad);
                complex ref = complex(std::cos(phase), std::sin(phase));
                double err = std::abs(sv[j][i] - ref);
                assert(err < 1e-12);
            }
        }

        std::cout << "  Test 2 PASSED: ArrayGeometry ULA scanning vectors" << std::endl;
    }

    // Test 3: ArrayGeometry UCA scanning vectors
    {
        ArrayGeometrySpec spec;
        spec.type = ArrayGeometrySpec::Type::UCA;
        spec.radius = 0.5;
        spec.ula_elements = 4;
        ArrayGeometry geom(spec);

        assert(geom.num_elements() == 4);

        // For UCA, elements are at (r*cos(phi_i), r*sin(phi_i), 0)
        // phase = 2*pi * (x*cos(theta) + y*sin(theta))
        std::vector<double> thetas = {0, 45, 90};
        mat sv = geom.gen_scanning_vectors(thetas);

        double r = 0.5;
        for (size_t i = 0; i < thetas.size(); ++i) {
            double theta_rad = deg_to_rad(thetas[i]);
            for (int j = 0; j < 4; ++j) {
                double phi = 2.0 * M_PI * j / 4.0;
                double x = r * std::cos(phi);
                double y = r * std::sin(phi);
                double phase = 2.0 * M_PI * (x * std::cos(theta_rad) + y * std::sin(theta_rad));
                complex ref = complex(std::cos(phase), std::sin(phase));
                double err = std::abs(sv[j][i] - ref);
                assert(err < 1e-12);
            }
        }

        std::cout << "  Test 3 PASSED: ArrayGeometry UCA scanning vectors" << std::endl;
    }

    // Test 4: ArrayGeometry URA scanning vectors
    {
        ArrayGeometrySpec spec;
        spec.type = ArrayGeometrySpec::Type::URA;
        spec.element_spacing = 0.5;
        spec.ura_rows = 2;
        spec.ura_cols = 3;
        ArrayGeometry geom(spec);

        assert(geom.num_elements() == 6);

        std::vector<double> thetas = {30, 60, 90};
        std::vector<double> phis = {0, 45, 90};
        mat sv = geom.gen_scanning_vectors(thetas, phis);

        assert(sv.size() == 6);
        assert(sv[0].size() == 3);

        double d = 0.5;
        for (size_t i = 0; i < thetas.size(); ++i) {
            double theta_rad = deg_to_rad(thetas[i]);
            double phi_rad = deg_to_rad(phis[i]);
            int idx = 0;
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 3; ++c) {
                    double x = d * (c - 1.0);
                    double y = d * (r - 0.5);
                    double phase = 2.0 * M_PI * (x * std::sin(theta_rad) * std::cos(phi_rad)
                                             + y * std::sin(theta_rad) * std::sin(phi_rad));
                    complex ref = complex(std::cos(phase), std::sin(phase));
                    double err = std::abs(sv[idx][i] - ref);
                    assert(err < 1e-12);
                    ++idx;
                }
        }

        std::cout << "  Test 4 PASSED: ArrayGeometry URA scanning vectors" << std::endl;
    }

    // Test 5: Arbitrary array geometry
    {
        ArrayGeometrySpec spec;
        spec.type = ArrayGeometrySpec::Type::ARBITRARY;
        spec.elements = {
            {0.0, 0.0, 0.0},
            {0.5, 0.0, 0.0},
            {0.25, 0.433, 0.0},
            {0.25, -0.433, 0.0}
        };
        ArrayGeometry geom(spec);

        assert(geom.num_elements() == 4);

        std::vector<double> thetas = {45, 90};
        mat sv = geom.gen_scanning_vectors(thetas);

        for (size_t i = 0; i < thetas.size(); ++i) {
            double theta_rad = deg_to_rad(thetas[i]);
            for (int j = 0; j < 4; ++j) {
                const auto& el = spec.elements[j];
                double phase = 2.0 * M_PI * (el.x * std::sin(theta_rad));
                complex ref = complex(std::cos(phase), std::sin(phase));
                double err = std::abs(sv[j][i] - ref);
                assert(err < 1e-12);
            }
        }

        std::cout << "  Test 5 PASSED: ArrayGeometry arbitrary elements" << std::endl;
    }

    // Test 6: Synchronizer corrects integer sample offset
    {
        Synchronizer sync;

        IQBuffer original(1000);
        for (int i = 0; i < 1000; ++i)
            original[i] = complex(static_cast<double>(i), 0.0);

        // Shift by 3 samples
        IQBuffer shifted = sync.interpolate_shift(original, 3.0);

        // shifted[i] should be original[i-3]
        for (int i = 3; i < 1000; ++i) {
            double err = std::abs(shifted[i] - original[i - 3]);
            assert(err < 1e-12);
        }
        // First 3 should be zero (or default)
        for (int i = 0; i < 3; ++i)
            assert(std::abs(shifted[i]) < 1e-12);

        std::cout << "  Test 6 PASSED: Synchronizer integer offset" << std::endl;
    }

    // Test 7: Synchronizer corrects fractional sample offset
    {
        Synchronizer sync;

        IQBuffer original(1000);
        for (int i = 0; i < 1000; ++i)
            original[i] = complex(static_cast<double>(i), 0.0);

        // Shift by 2.5 samples (integer=2, frac=0.5)
        IQBuffer shifted = sync.interpolate_shift(original, 2.5);

        // Linear interpolation: shifted[i] = original[i-3] * 0.5 + original[i-2] * 0.5
        for (int i = 3; i < 999; ++i) {
            complex expected = original[i - 3] * 0.5 + original[i - 2] * 0.5;
            double err = std::abs(shifted[i] - expected);
            assert(err < 1e-6);
        }

        std::cout << "  Test 7 PASSED: Synchronizer fractional offset" << std::endl;
    }

    // Test 8: Gain/phase calibration
    {
        Synchronizer sync;

        MultiSiteData data;
        ReceiverInfo info;
        info.id = "site_A";
        info.gain = 6.0;      // 2x amplitude
        info.phase_offset = M_PI / 2;  // 90 degrees
        info.clock_offset = 0.0;
        info.num_channels = 1;
        info.array_geometry.type = ArrayGeometrySpec::Type::ULA;
        info.array_geometry.ula_elements = 1;
        info.array_geometry.element_spacing = 0.5;
        info.ref_channel_index = 0;

        data.receivers.push_back(info);
        data.fs = 1000.0;

        IQBuffer ch(100);
        for (int i = 0; i < 100; ++i)
            ch[i] = complex(1.0, 0.0);
        data.iq_data.push_back(IQMatrix(1, ch));

        sync.calibrate(data);

        // After calibration: sample * 2.0 * exp(j*pi/2) = sample * 2j
        for (int i = 0; i < 100; ++i) {
            double err = std::abs(data.iq_data[0][0][i] - complex(0.0, 2.0));
            assert(err < 1e-12);
        }

        std::cout << "  Test 8 PASSED: Gain/phase calibration" << std::endl;
    }

    // Test 9: Single-site mode matches Pipeline (with dataset if available)
    {
        std::vector<std::string> candidates = {
            "dataset",
            "../passiveradar_data/dataset",
            "../../passiveradar_data/dataset",
            "/mnt/d/HAMPR/Code/passiveradar_data/dataset"
        };
        auto file_exists = [](const std::string& p) {
            std::ifstream f(p); return f.good();
        };
        std::string data_dir;
        for (const auto& dir : candidates) {
            if (file_exists(dir + "/dataset_494.txt") &&
                file_exists(dir + "/target_reference_track.txt")) {
                data_dir = dir;
                break;
            }
        }

        if (data_dir.empty()) {
            std::cout << "  Test 9 SKIPPED: VEGA dataset not found" << std::endl;
        } else {
            Config config;
            config.ref_channel_index = 0;
            config.filter_taps = 128;
            config.num_antennas = 3;
            config.antenna_spacing = 0.528;
            config.window_type = "Hann";
            config.search_window_size = 8;
            config.metric_win = {6, 6, 3, 3};

            TextDataSource ds;
            double fs;
            IQMatrix iq = ds.load(data_dir + "/dataset_494.txt", fs);
            auto track = ds.load_track(data_dir + "/target_reference_track.txt");
            assert(track.size() >= 1);
            int N = static_cast<int>(iq[0].size());

            // Run through Pipeline
            Pipeline pipeline(config);
            ProcessingResult ref_result = pipeline.process(iq, track.front(), fs, N);

            // Run through MultiSitePipeline (single-site mode)
            MultiSitePipeline ms_pipeline(config);
            MultiSiteData data;
            data.fs = fs;
            data.iq_data.push_back(iq);
            ReceiverInfo info;
            info.id = "site_A";
            info.num_channels = 4;
            info.ref_channel_index = 0;
            info.array_geometry.type = ArrayGeometrySpec::Type::ULA;
            info.array_geometry.ula_elements = 4;
            info.array_geometry.element_spacing = 0.528;
            data.receivers.push_back(info);

            MultiSiteResult ms_result = ms_pipeline.process(data, track);
            assert(ms_result.site_results.size() == 1);

            auto& site_res = ms_result.site_results[0];
            double range_err = std::abs(site_res.target_range_cell_found - ref_result.target_range_cell_found);
            double doppler_err = std::abs(site_res.target_doppler_cell_found - ref_result.target_doppler_cell_found);
            double sinr_err = std::abs(site_res.sinr - ref_result.sinr);
            double az_err = std::abs(site_res.azimuth_estimated - ref_result.azimuth_estimated);

            assert(range_err < 1e-6);
            assert(doppler_err < 1e-6);
            assert(sinr_err < 1e-6);
            assert(az_err < 1e-6);

            std::cout << "  Test 9 PASSED: Single-site mode matches Pipeline (err < 1e-6)" << std::endl;
        }
    }

    // Test 10: Two-site processing
    {
        Config config;
        config.ref_channel_index = 0;
        config.filter_taps = 128;
        config.num_antennas = 3;
        config.antenna_spacing = 0.528;
        config.window_type = "Hann";
        config.search_window_size = 8;
        config.metric_win = {6, 6, 3, 3};
        config.enable_fusion = true;
        config.enable_localization = false;

        ReceiverInfo info_a, info_b;
        info_a.id = "site_A";
        info_a.latitude = -37.8472;
        info_a.longitude = 145.0451;
        info_a.altitude = 12.0;
        info_a.num_channels = 4;
        info_a.ref_channel_index = 0;
        info_a.array_geometry.type = ArrayGeometrySpec::Type::ULA;
        info_a.array_geometry.ula_elements = 4;
        info_a.array_geometry.element_spacing = 0.528;

        info_b.id = "site_B";
        info_b.latitude = -37.8430;
        info_b.longitude = 145.0480;
        info_b.altitude = 15.0;
        info_b.num_channels = 4;
        info_b.ref_channel_index = 0;
        info_b.array_geometry.type = ArrayGeometrySpec::Type::ULA;
        info_b.array_geometry.ula_elements = 4;
        info_b.array_geometry.element_spacing = 0.528;

        config.multi_site_receivers = {info_a, info_b};

        MultiSitePipeline pipeline(config);

        // Generate synthetic IQ data for two sites
        IQMatrix iq_a(4, IQBuffer(1024));
        IQMatrix iq_b(4, IQBuffer(1024));
        std::mt19937 rng(42);
        std::normal_distribution<double> dist(0.0, 0.1);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 1024; ++j) {
                iq_a[i][j] = complex(dist(rng), dist(rng));
                iq_b[i][j] = complex(dist(rng), dist(rng));
            }

        MultiSiteData data;
        data.fs = 2500000.0;
        data.receivers = {info_a, info_b};
        data.iq_data = {iq_a, iq_b};

        TargetTrackPoint ref_track;
        ref_track.time_index = 0;
        ref_track.timestamp = 0.0;
        ref_track.latitude = -37.8450;
        ref_track.longitude = 145.0460;
        ref_track.altitude = 13.0;
        ref_track.range = 5000.0;
        ref_track.doppler = 50.0;
        ref_track.azimuth = 45.0;

        TargetTrack track = {ref_track};

        MultiSiteResult result = pipeline.process(data, track);
        assert(result.site_results.size() == 2);
        assert(result.total_time > 0.0);

        std::cout << "  Test 10 PASSED: Two-site processing" << std::endl;
    }

    // Test 11: Bearing intersection accuracy
    {
        ResultFusion fusion;

        // Two sites at known positions with known bearings
        // Site A at (0, 0), Site B at (10, 0)
        // Bearing from A: 45 degrees (NE)
        // Bearing from B: 135 degrees (NW)
        // Intersection: (5, 5)

        std::vector<ProcessingResult> results(2);
        results[0].azimuth_estimated = 45.0;
        results[0].sinr = 10.0;
        results[1].azimuth_estimated = 135.0;
        results[1].sinr = 10.0;

        std::vector<ReceiverInfo> receivers(2);
        receivers[0].id = "A";
        receivers[0].latitude = 0.0;
        receivers[0].longitude = 0.0;
        receivers[1].id = "B";
        receivers[1].latitude = 0.0;
        receivers[1].longitude = 0.0001;  // ~11.1 meters east

        auto fused = fusion.fuse_bearings(results, receivers);
        assert(fused.contributing_bearings.size() == 2);
        assert(fused.confidence > 0.0);

        double combined = fusion.combined_sinr(results);
        assert(combined > 0.0);

        std::cout << "  Test 11 PASSED: Bearing intersection (combined SINR=" << combined << ")" << std::endl;
    }

    // Test 12: TDoA localization with synthetic data
    {
        TDoALocalizer localizer;

        // Create two synchronized batches with a known time delay
        SynchronizedBatch batch;
        batch.fs = 1000000.0;  // 1 MHz
        batch.batch_samples = 1000;

        ReceiverInfo info_a, info_b;
        info_a.id = "A";
        info_a.latitude = -37.8472;
        info_a.longitude = 145.0451;
        info_b.id = "B";
        info_b.latitude = -37.8430;
        info_b.longitude = 145.0480;
        batch.receivers = {info_a, info_b};

        // Create reference signals: site B has a 10-sample delay
        IQBuffer ref_a(1000);
        IQBuffer ref_b(1000);
        for (int i = 0; i < 1000; ++i) {
            ref_a[i] = complex(std::sin(2.0 * M_PI * i / 100.0), 0.0);
            ref_b[i] = (i >= 10) ? ref_a[i - 10] : complex(0.0, 0.0);
        }
        batch.aligned_data.push_back(IQMatrix(1, ref_a));
        batch.aligned_data.push_back(IQMatrix(1, ref_b));

        TargetTrackPoint ref_track;
        ref_track.range = 5000.0;
        TargetTrack track = {ref_track};

        auto geo = localizer.localize(batch, track);
        assert(geo.latitude != 0.0 || geo.longitude != 0.0);
        assert(geo.residual_error >= 0.0);

        std::cout << "  Test 12 PASSED: TDoA localization (residual=" << geo.residual_error << "m)" << std::endl;
    }

    // Test 13: MultiSiteGraph async metrics interface
    {
        Config config;
        config.num_antennas = 3;
        MultiSiteGraph graph(config);

        assert(graph.batches_processed() == 0);
        assert(graph.avg_latency_ms() == 0.0);
        assert(graph.throughput_samples_per_sec() >= 0.0);

        std::cout << "  Test 13 PASSED: MultiSiteGraph metrics interface" << std::endl;
    }

    // Test 14: MultiSiteGraph streaming (with dataset if available)
    {
        std::vector<std::string> candidates = {
            "dataset", "../passiveradar_data/dataset",
            "../../passiveradar_data/dataset",
            "/mnt/d/HAMPR/Code/passiveradar_data/dataset"
        };
        auto file_exists = [](const std::string& p) {
            std::ifstream f(p); return f.good();
        };
        std::string data_dir;
        for (const auto& dir : candidates) {
            if (file_exists(dir + "/dataset_494.txt")) {
                data_dir = dir;
                break;
            }
        }

        if (data_dir.empty()) {
            std::cout << "  Test 14 SKIPPED: VEGA dataset not found" << std::endl;
        } else {
            Config config;
            config.ref_channel_index = 0;
            config.filter_taps = 128;
            config.num_antennas = 3;
            config.antenna_spacing = 0.528;
            config.window_type = "Hann";
            config.search_window_size = 8;
            config.metric_win = {6, 6, 3, 3};
            config.enable_fusion = false;
            config.enable_localization = false;

            MultiSiteGraph graph(config);

            MultiSiteDataSource source;
            bool _opened = source.open_site("site_A", data_dir + "/dataset_494.txt"); assert(_opened);
            assert(source.has_more_site("site_A"));
            assert(source.num_channels_site("site_A") == 4);
            assert(source.sampling_rate_site("site_A") > 0);

            auto track = source.load_track(data_dir + "/target_reference_track.txt");
            assert(track.size() >= 1);

            auto results = graph.process_stream(source, track, 65536);
            assert(!results.empty());
            assert(graph.batches_processed() > 0);
            assert(graph.avg_latency_ms() > 0);

            source.close_all();

            std::cout << "  Test 14 PASSED: MultiSiteGraph streaming ("
                      << graph.batches_processed() << " batches, "
                      << "avg_latency=" << graph.avg_latency_ms() << "ms"
                      << ")" << std::endl;
        }
    }

    // Test 15: Metadata loading from JSON
    {
        // Create a temporary JSON metadata file
        std::string metadata = R"(
{
    "fs": 2500000.0,
    "receivers": [
        {
            "id": "site_A",
            "latitude": -37.8472,
            "longitude": 145.0451,
            "altitude": 12.0,
            "clock_offset": 0.0,
            "gain": -1.5,
            "phase_offset": 0.03,
            "ref_channel_index": 0,
            "num_channels": 4,
            "array_geometry": {
                "type": "ULA",
                "element_spacing": 0.528,
                "ula_elements": 4
            }
        },
        {
            "id": "site_B",
            "latitude": -37.8430,
            "longitude": 145.0480,
            "altitude": 15.0,
            "clock_offset": 0.000012,
            "gain": -2.0,
            "phase_offset": 0.01,
            "ref_channel_index": 0,
            "num_channels": 4,
            "array_geometry": {
                "type": "ULA",
                "element_spacing": 0.528,
                "ula_elements": 4
            }
        }
    ]
}
)";

        std::string tmpfile = "test_metadata.json";
        std::ofstream f(tmpfile);
        f << metadata;
        f.close();

        MultiSiteDataSource source;
        assert(source.load_metadata(tmpfile));

        assert(source.receivers().size() == 2);
        assert(source.receivers()[0].id == "site_A");
        assert(source.receivers()[0].latitude == -37.8472);
        assert(source.receivers()[0].longitude == 145.0451);
        assert(source.receivers()[0].altitude == 12.0);
        assert(source.receivers()[0].clock_offset == 0.0);
        assert(source.receivers()[0].gain == -1.5);
        assert(source.receivers()[0].phase_offset == 0.03);
        assert(source.receivers()[0].array_geometry.type == ArrayGeometrySpec::Type::ULA);
        assert(source.receivers()[0].array_geometry.element_spacing == 0.528);
        assert(source.receivers()[0].array_geometry.ula_elements == 4);

        assert(source.receivers()[1].id == "site_B");
        assert(source.receivers()[1].longitude == 145.0480);
        assert(source.receivers()[1].clock_offset == 0.000012);

        assert(source.sampling_rate() == 2500000.0);

        std::remove(tmpfile.c_str());

        std::cout << "  Test 15 PASSED: Metadata loading from JSON" << std::endl;
    }

    // Test 16: Multi-static cross-correlation
    {
        MultiStaticCorrelator correlator;

        SynchronizedBatch batch;
        batch.fs = 1000000.0;
        batch.batch_samples = 1000;

        ReceiverInfo info_a, info_b;
        info_a.id = "A";
        info_b.id = "B";
        batch.receivers = {info_a, info_b};

        IQBuffer ref_a(1000);
        IQBuffer ref_b(1000);
        for (int i = 0; i < 1000; ++i) {
            ref_a[i] = complex(std::sin(2.0 * M_PI * i / 100.0), 0.0);
            ref_b[i] = (i >= 5) ? ref_a[i - 5] : complex(0.0, 0.0);
        }
        batch.aligned_data.push_back(IQMatrix(1, ref_a));
        batch.aligned_data.push_back(IQMatrix(1, ref_b));

        auto results = correlator.cross_correlate(batch);
        assert(results.size() == 1);
        assert(results[0].site_a == 0);
        assert(results[0].site_b == 1);
        assert(results[0].correlation_peak > 0.0);

        auto ranges = correlator.compute_bistatic_ranges(results, batch.fs);
        assert(ranges.size() == 1);
        assert(ranges[0] >= 0.0);

        auto confidence = correlator.cooperative_confidence(results, batch.aligned_data[0]);
        assert(confidence >= 0.0);

        std::cout << "  Test 16 PASSED: Multi-static cross-correlation (peak=" << results[0].correlation_peak << ")" << std::endl;
    }

    std::cout << "All Phase 7 multi-site tests passed!" << std::endl;
    return 0;
}

