// numerical_dump.cpp - Stage-by-stage intermediate dump for HAMPR.
//
// Mirrors the exact stage wiring of hampr::Pipeline::process (see
// src/dsp/pipeline.cpp) and writes each DSP stage's intermediate output to a
// portable binary format so it can be compared against the benchmark
// reference dump (validation/bench_dump.cpp) by tools/compare_dumps.
//
// Dump format (little-endian):
//   bytes  0..3 : magic "HDMP"
//   int32       : rows
//   int32       : cols
//   rows*cols   : complex<double> (real, imag)   -- written row by row
#include <hampr/dsp/channel_isolator.hpp>
#include <hampr/dsp/beamformer.hpp>
#include <hampr/dsp/clutter_canceler.hpp>
#include <hampr/dsp/detector.hpp>
#include <hampr/dsp/doa_estimator.hpp>
#include <hampr/dsp/metric_extractor.hpp>
#include <hampr/io/text_data_source.hpp>
#include <hampr/core/config.hpp>
#include <hampr/utils/math_utils.hpp>

#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace hampr;

static bool file_exists(const std::string& p) {
    std::ifstream f(p);
    return f.good();
}

static std::string find_dataset_dir() {
    const std::vector<std::string> cand = {
        "dataset",
        "../passiveradar_data/dataset",
        "../../passiveradar_data/dataset",
        "Code/passiveradar_data/dataset",
        "/mnt/d/HAMPR/Code/passiveradar_data/dataset",
    };
    for (const auto& d : cand) {
        if (file_exists(d + "/target_reference_track.txt") &&
            file_exists(d + "/dataset_494.txt"))
            return d;
    }
    return "dataset";
}

static void write_header(std::ofstream& f, int32_t rows, int32_t cols) {
    const char magic[4] = {'H', 'D', 'M', 'P'};
    f.write(magic, 4);
    f.write(reinterpret_cast<const char*>(&rows), 4);
    f.write(reinterpret_cast<const char*>(&cols), 4);
}

static void dump_1d(const std::string& path, const IQBuffer& a) {
    std::ofstream f(path, std::ios::binary);
    write_header(f, 1, static_cast<int32_t>(a.size()));
    f.write(reinterpret_cast<const char*>(a.data()), a.size() * sizeof(complex));
}

static void dump_2d(const std::string& path, const mat& m) {
    std::ofstream f(path, std::ios::binary);
    const int32_t rows = static_cast<int32_t>(m.size());
    const int32_t cols = m.empty() ? 0 : static_cast<int32_t>(m[0].size());
    write_header(f, rows, cols);
    for (const auto& row : m)
        f.write(reinterpret_cast<const char*>(row.data()), cols * sizeof(complex));
}

static void dump_weights(const std::string& path, const std::vector<complex>& w) {
    std::ofstream f(path, std::ios::binary);
    write_header(f, static_cast<int32_t>(w.size()), 1);
    f.write(reinterpret_cast<const char*>(w.data()), w.size() * sizeof(complex));
}

static void dump_results(const std::string& path,
                         double ref_range, double range_found,
                         double ref_doppler, double doppler_found,
                         double az_ref, double az_est, double sinr) {
    const complex c[7] = {
        complex(ref_range, 0.0), complex(range_found, 0.0),
        complex(ref_doppler, 0.0), complex(doppler_found, 0.0),
        complex(az_ref, 0.0), complex(az_est, 0.0), complex(sinr, 0.0)};
    std::ofstream f(path, std::ios::binary);
    write_header(f, 1, 7);
    f.write(reinterpret_cast<const char*>(c), 7 * sizeof(complex));
}

int main(int argc, char** argv) {
    const std::string out_root = (argc > 1) ? argv[1] : "../val_results";
    const std::string data_dir = find_dataset_dir();
    std::cerr << "[hampr_dump] data_dir=" << data_dir
              << " out_root=" << out_root << std::endl;

    Config config;
    config.ref_channel_index = 0;
    config.filter_taps = 128;
    config.num_antennas = 3;
    config.antenna_spacing = 0.528;
    config.window_type = "Hann";
    config.search_window_size = 8;
    config.metric_win = {6, 6, 3, 3};

    ChannelIsolator isolator(config.ref_channel_index);
    Beamformer beamformer(config.antenna_spacing);
    ClutterCanceler canceler(config.filter_taps);
    Detector detector(config.window_type);
    DOAEstimator doa_estimator;
    MetricExtractor metric_extractor;

    TextDataSource data_source;
    const TargetTrack track = data_source.load_track(data_dir + "/target_reference_track.txt");

    for (const auto& ref_track : track) {
        const int idx = ref_track.time_index;
        const std::string fname = data_dir + "/dataset_" + std::to_string(idx) + ".txt";
        double fs = 0.0;
        const IQMatrix iq_data = data_source.load(fname, fs);
        const int N = static_cast<int>(iq_data[0].size());
        const std::string dir = out_root + "/dump/" + std::to_string(idx);
        std::filesystem::create_directories(dir);

        const double max_range = ref_track.range;
        const double max_doppler = ref_track.doppler;
        std::vector<double> alignment(config.num_antennas);
        for (int i = 0; i < config.num_antennas; ++i)
            alignment[i] = i * config.antenna_spacing;
        const double fD_res = fs / (2.0 * N);
        const double range_res = 3.0e8 / fs;
        const int doppler_cells =
            static_cast<int>(max_doppler / fD_res) + config.search_window_size +
            config.metric_win[0] + 1;
        const int max_doppler_ext = static_cast<int>(std::ceil(doppler_cells * fD_res));
        const int max_range_ext = static_cast<int>(
            std::pow(2, std::ceil(std::log2(max_range / range_res +
                                            config.search_window_size + config.metric_win[0]))));
        const int ref_range_cell = static_cast<int>(std::round(ref_track.range / range_res));
        const int ref_doppler_cell =
            static_cast<int>(std::round(ref_track.doppler / fD_res)) + doppler_cells;
        double doa_dir = ref_track.azimuth;
        if (doa_dir < 0) doa_dir += 180;
        if (doa_dir > 180) doa_dir -= 180;

        // Stage 1: Channel isolation
        auto ref_and_survs = isolator.isolate(iq_data);
        IQBuffer& ref_ch = ref_and_survs.first;
        const IQMatrix& surv_chs = ref_and_survs.second;
        dump_1d(dir + "/01_ref_ch.bin", ref_ch);
        for (size_t m = 0; m < surv_chs.size(); ++m)
            dump_1d(dir + "/02_surv_ch_" + std::to_string(m) + ".bin", surv_chs[m]);

        // Stage 2: Beamforming
        const auto surv_chs_t = transpose(surv_chs);
        const auto beam_result = beamformer.beamform(surv_chs_t, alignment, doa_dir);
        const IQBuffer surv_ch = beam_result.surveillance_channel;
        dump_weights(dir + "/03_beam_w.bin", beam_result.weights);
        dump_1d(dir + "/03_beam_surv.bin", surv_ch);

        // Stage 3: Clutter cancellation (Wiener SMI-MRE)
        const IQBuffer filt_surv_ch = canceler.filter(ref_ch, surv_ch);
        dump_1d(dir + "/04_clutter.bin", filt_surv_ch);

        // Stage 4: Windowing
        IQBuffer windowed = detector.apply_window(filt_surv_ch);
        dump_1d(dir + "/05_windowed.bin", windowed);

        // Stage 5: Detection (cross-correlation range-Doppler map)
        RDMap rd_matrix = detector.detect(ref_ch, windowed, fs, max_doppler_ext, max_range_ext);
        dump_2d(dir + "/06_rdmap.bin", rd_matrix);

        // Stage 6: Target finding (peak search in window)
        int range_found = 0, doppler_found = 0;
        double max_val = 0.0;
        for (auto& row : rd_matrix)
            for (auto& c : row)
                c = complex(std::abs(c), 0.0);
        for (int i = ref_doppler_cell - config.search_window_size;
             i < ref_doppler_cell + config.search_window_size; ++i)
            for (int j = ref_range_cell - config.search_window_size;
                 j < ref_range_cell + config.search_window_size; ++j)
                if (rd_matrix[i][j].real() > max_val) {
                    max_val = rd_matrix[i][j].real();
                    doppler_found = i;
                    range_found = j;
                }

        // Stage 7: DOA - per-channel RD maps + MUSIC azimuth
        const std::vector<std::pair<int, int>> hit_list = {{range_found, doppler_found}};
        std::vector<RDMap> rd_maps_for_doa(surv_chs.size());
        for (size_t m = 0; m < surv_chs.size(); ++m) {
            IQBuffer f = canceler.filter(ref_ch, surv_chs[m]);
            IQBuffer w = detector.apply_window(f);
            rd_maps_for_doa[m] = detector.detect(ref_ch, w, fs, max_doppler_ext, max_range_ext);
            dump_2d(dir + "/07_doa_rd_" + std::to_string(m) + ".bin", rd_maps_for_doa[m]);
        }
        const auto doa = doa_estimator.estimate(rd_maps_for_doa, hit_list, alignment);
        const double azimuth_est = doa.empty() ? 0.0 : doa[0];

        // Stage 8: Metric extraction (SNR)
        std::vector<int> win = config.metric_win;
        std::vector<int> target_rd = {range_found, doppler_found};
        target_rd[1] -= (static_cast<int>(rd_matrix.size()) - 1) / 2;
        const double sinr = metric_extractor.extract_snr(rd_matrix, target_rd, win);

        dump_results(dir + "/08_results.bin",
                     static_cast<double>(ref_range_cell), static_cast<double>(range_found),
                     static_cast<double>(ref_doppler_cell), static_cast<double>(doppler_found),
                     ref_track.azimuth, azimuth_est, sinr);

        std::cout << "idx=" << idx << " range " << ref_range_cell << "/" << range_found
                  << " doppler " << ref_doppler_cell << "/" << doppler_found
                  << " az " << ref_track.azimuth << "/" << azimuth_est
                  << " sinr " << sinr << std::endl;
    }
    return 0;
}
