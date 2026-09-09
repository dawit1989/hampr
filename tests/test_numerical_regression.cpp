#include <hampr/dsp/pipeline.hpp>
#include <hampr/io/text_data_source.hpp>
#include <hampr/core/config.hpp>
#include <hampr/core/exception.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <iomanip>

static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

static std::string find_dataset_dir() {
    std::vector<std::string> candidates = {
        "dataset",
        "../passiveradar_data/dataset",
        "../../passiveradar_data/dataset",
        "/mnt/d/HAMPR/Code/passiveradar_data/dataset"
    };
    for (const auto& dir : candidates) {
        if (file_exists(dir + "/target_reference_track.txt") &&
            file_exists(dir + "/dataset_494.txt")) {
            return dir;
        }
    }
    return "dataset";
}

struct ExpectedResult {
    int range_cell_ref;
    int range_cell_found;
    int doppler_cell_ref;
    int doppler_cell_found;
    double azimuth_ref;
    double azimuth_estimated;
    double sinr;
};

static std::vector<ExpectedResult> get_reference_results() {
    return {
        {103, 110, 168, 167, 84.4083, 72, 9.1824},
        {103,  99, 168, 170, 84.4083, 67, 7.6557},
        {103, 109, 168, 160, 84.4083, 67, 9.00417},
        {103,  98, 168, 172, 84.4083, 66, 7.05235},
        {103, 106, 168, 172, 84.4083, 66, 9.02792}
    };
}

int main() {
    using namespace hampr;

    std::string data_dir = find_dataset_dir();
    std::cout << "Numerical regression test: using dataset from " << data_dir << std::endl;

    Config config;
    config.ref_channel_index = 0;
    config.filter_taps = 128;
    config.num_antennas = 3;
    config.antenna_spacing = 0.528;
    config.window_type = "Hann";
    config.search_window_size = 8;
    config.metric_win = {6, 6, 3, 3};

    Pipeline pipeline(config);
    TextDataSource data_source;

    auto track = data_source.load_track(data_dir + "/target_reference_track.txt");
    assert(track.size() >= 5);

    auto expected = get_reference_results();

    int test_idx = 0;
    double max_snr_err = 0.0;
    double max_range_err = 0.0;
    double max_doppler_err = 0.0;
    double max_az_err = 0.0;
    bool all_passed = true;

    for (size_t t = 0; t < track.size() && test_idx < 5; ++t) {
        int idx = track[t].time_index;
        std::string filename = data_dir + "/dataset_" + std::to_string(idx) + ".txt";
        if (!file_exists(filename)) continue;

        double fs;
        auto iq_data = data_source.load(filename, fs);
        int N = static_cast<int>(iq_data[0].size());

        auto result = pipeline.process(iq_data, track[t], fs, N);

        const auto& exp = expected[test_idx];

        int range_err = std::abs(result.target_range_cell_found - exp.range_cell_found);
        int doppler_err = std::abs(result.target_doppler_cell_found - exp.doppler_cell_found);
        double snr_err = std::abs(result.sinr - exp.sinr);
        double az_err = std::abs(result.azimuth_estimated - exp.azimuth_estimated);
        if (az_err > 180.0) az_err = 360.0 - az_err;

        max_snr_err = std::max(max_snr_err, snr_err);
        max_range_err = std::max(max_range_err, static_cast<double>(range_err));
        max_doppler_err = std::max(max_doppler_err, static_cast<double>(doppler_err));
        max_az_err = std::max(max_az_err, az_err);

        bool stage_pass = (range_err <= 10) && (doppler_err <= 20) &&
                          (snr_err < 5.0) && (az_err <= 90.0);

        if (!stage_pass) all_passed = false;

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Dataset " << idx << ":" << std::endl;
        std::cout << "  Range cell:   ref=" << exp.range_cell_found
                  << "  hampr=" << result.target_range_cell_found
                  << "  diff=" << range_err << std::endl;
        std::cout << "  Doppler cell: ref=" << exp.doppler_cell_found
                  << "  hampr=" << result.target_doppler_cell_found
                  << "  diff=" << doppler_err << std::endl;
        std::cout << "  Azimuth:      ref=" << exp.azimuth_estimated
                  << "  hampr=" << result.azimuth_estimated
                  << "  diff=" << az_err << std::endl;
        std::cout << "  SINR:         ref=" << exp.sinr
                  << "  hampr=" << result.sinr
                  << "  diff=" << snr_err << std::endl;

        test_idx++;
    }

    std::cout << std::endl;
    std::cout << "=== Summary ===" << std::endl;
    std::cout << "Datasets tested: " << test_idx << std::endl;
    std::cout << "Max range cell error:    " << max_range_err << std::endl;
    std::cout << "Max doppler cell error:  " << max_doppler_err << std::endl;
    std::cout << "Max azimuth error:       " << max_az_err << std::endl;
    std::cout << "Max SINR error:          " << max_snr_err << std::endl;

    if (all_passed) {
        std::cout << "All numerical regression assertions passed!" << std::endl;
        return 0;
    } else {
        std::cout << "Some numerical regression assertions FAILED!" << std::endl;
        return 1;
    }
}