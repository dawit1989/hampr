#include <hampr/dsp/pipeline.hpp>
#include <hampr/io/text_data_source.hpp>
#include <hampr/core/config.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>

static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

static std::string find_dataset_dir() {
    std::vector<std::string> candidates = {
        "dataset",
        "../passiveradar_data/dataset",
        "../../passiveradar_data/dataset",
        "/mnt/d/HAMPR/Code/passiveradar_data/dataset",
    };
    for (const auto& dir : candidates) {
        if (file_exists(dir + "/target_reference_track.txt") &&
            file_exists(dir + "/dataset_494.txt")) {
            return dir;
        }
    }
    return "dataset";
}

int main() {
    using namespace hampr;

    std::string data_dir = find_dataset_dir();
    std::cout << "Regression test: using dataset from " << data_dir << std::endl;

    hampr::Config config;
    config.ref_channel_index = 0;
    config.filter_taps = 128;
    config.num_antennas = 3;
    config.antenna_spacing = 0.528;
    config.window_type = "Hann";
    config.search_window_size = 8;
    config.metric_win = {6, 6, 3, 3};

    hampr::Pipeline pipeline(config);
    hampr::TextDataSource data_source;

    auto track = data_source.load_track(data_dir + "/target_reference_track.txt");
    assert(track.size() >= 1);

    int idx = track.front().time_index;
    double fs;
    std::string filename = data_dir + "/dataset_" + std::to_string(idx) + ".txt";
    auto iq_data = data_source.load(filename, fs);
    int N = static_cast<int>(iq_data[0].size());

    auto result = pipeline.process(iq_data, track.front(), fs, N);

    std::cout << "Results for index " << idx << ":" << std::endl;
    std::cout << "  range cell: ref=" << result.target_range_cell_ref
              << " found=" << result.target_range_cell_found << std::endl;
    std::cout << "  doppler cell: ref=" << result.target_doppler_cell_ref
              << " found=" << result.target_doppler_cell_found << std::endl;
    std::cout << "  azimuth: ref=" << result.azimuth_ref
              << " estimated=" << result.azimuth_estimated << std::endl;
    std::cout << "  SINR: " << result.sinr << std::endl;
    std::cout << "  total_time: " << result.total_time << "s" << std::endl;

    int range_diff = std::abs(result.target_range_cell_found - result.target_range_cell_ref);
    assert(range_diff <= 10);

    int doppler_diff = std::abs(result.target_doppler_cell_found - result.target_doppler_cell_ref);
    assert(doppler_diff <= 20);

    assert(result.sinr > 0.0 && result.sinr < 30.0);

    double az_diff = std::abs(result.azimuth_estimated - result.azimuth_ref);
    if (az_diff > 180.0) az_diff = 360.0 - az_diff;
    assert(az_diff <= 90.0);

    std::cout << "All regression assertions passed!" << std::endl;
    return 0;
}

