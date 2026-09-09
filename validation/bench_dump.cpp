// bench_dump.cpp - Stage-by-stage intermediate dump for the benchmark.
//
// This is the AUTHORITY reference driver. It links the unmodified benchmark
// object files (channelPreparation.o, beamform.o, clutterCancellation.o,
// detector.o, directionEstimation.o, metricExtract.o, helper.o) and exercises
// each DSP stage function exactly as the benchmark main.cpp does, writing the
// same portable binary dump format as tests/numerical_dump.cpp so the two
// implementations can be compared stage-by-stage by tools/compare_dumps.
//
// It does NOT modify the benchmark source; it only adds a separate harness.
#include <algorithm>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "channelPreparation.h"
#include "beamform.h"
#include "clutterCancellation.h"
#include "detector.h"
#include "directionEstimation.h"
#include "metricExtract.h"
#include "helper.h"
#include "numpy_helper.h"

typedef std::vector<double> double_array;
typedef std::vector<double_array> double_mat;

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

static void load_target_track(double_mat& track, const std::string& path) {
    std::ifstream infile(path);
    int rows, cols;
    infile >> rows >> cols;
    track = double_mat(rows, double_array(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            infile >> track[i][j];
}

static void load_dataset(const std::string& filename, mat& data, double& fs) {
    std::ifstream infile(filename);
    infile >> fs;
    int rows, cols;
    infile >> rows >> cols;
    data = mat(rows, array(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) {
            double re, im;
            infile >> re >> im;
            data[i][j] = std::complex<double>(re, im);
        }
}

static void write_header(std::ofstream& f, int32_t rows, int32_t cols) {
    const char magic[4] = {'H', 'D', 'M', 'P'};
    f.write(magic, 4);
    f.write(reinterpret_cast<const char*>(&rows), 4);
    f.write(reinterpret_cast<const char*>(&cols), 4);
}

static void dump_1d(const std::string& path, const array& a) {
    std::ofstream f(path, std::ios::binary);
    write_header(f, 1, static_cast<int32_t>(a.size()));
    f.write(reinterpret_cast<const char*>(a.data()), a.size() * sizeof(std::complex<double>));
}

static void dump_2d(const std::string& path, const mat& m) {
    std::ofstream f(path, std::ios::binary);
    const int32_t rows = static_cast<int32_t>(m.size());
    const int32_t cols = m.empty() ? 0 : static_cast<int32_t>(m[0].size());
    write_header(f, rows, cols);
    for (const auto& row : m)
        f.write(reinterpret_cast<const char*>(row.data()), cols * sizeof(std::complex<double>));
}

static void dump_weights(const std::string& path, const mat& w) {
    // w is M x 1
    std::ofstream f(path, std::ios::binary);
    write_header(f, static_cast<int32_t>(w.size()), 1);
    for (const auto& row : w)
        f.write(reinterpret_cast<const char*>(row.data()), sizeof(std::complex<double>));
}

static void dump_results(const std::string& path,
                         double ref_range, double range_found,
                         double ref_doppler, double doppler_found,
                         double az_ref, double az_est, double sinr) {
    const std::complex<double> c[7] = {
        {ref_range, 0.0}, {range_found, 0.0}, {ref_doppler, 0.0},
        {doppler_found, 0.0}, {az_ref, 0.0}, {az_est, 0.0}, {sinr, 0.0}};
    std::ofstream f(path, std::ios::binary);
    write_header(f, 1, 7);
    f.write(reinterpret_cast<const char*>(c), 7 * sizeof(std::complex<double>));
}

static mat estimate_corr_matrix(const mat& X) {
    const std::complex<double> N(static_cast<double>(X.size()), 0.0);
    const int M = static_cast<int>(X[0].size());
    mat R(M, array(M));
    mat t = transpose(X);
    R = matMult(t, transpose(conjugate(t)));
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            R[i][j] /= N;
    return R;
}

int main(int argc, char** argv) {
    const std::string out_root = (argc > 1) ? argv[1] : "../val_results";
    const std::string data_dir = find_dataset_dir();
    std::cerr << "[bench_dump] data_dir=" << data_dir
              << " out_root=" << out_root << std::endl;

    const int ref_ch_index = 0;
    const int K = 128;
    const int M = 3;
    const double d = 0.528;
    const std::string cc_det_window = "Hann";
    std::vector<int> win = {6, 6, 3, 3};
    std::vector<int> win_pos = {30, 50};
    const int search_window_size = 8;

    double_mat target_reference_track;
    load_target_track(target_reference_track, data_dir + "/target_reference_track.txt");

    const int start_file_idx = static_cast<int>(target_reference_track.front()[0]);
    const int stop_file_idx = static_cast<int>(target_reference_track.back()[0]);

    double max_range = -std::numeric_limits<double>::infinity();
    double max_Doppler = max_range;
    for (const auto& row : target_reference_track) {
        max_range = std::max(max_range, row[7]);
        max_Doppler = std::max(max_Doppler, row[8]);
    }

    std::vector<double> array_alignment(M);
    for (int i = 0; i < M; ++i)
        array_alignment[i] = i * d;

    for (int t = 0; t < stop_file_idx - start_file_idx + 1; ++t) {
        const int idx = t + start_file_idx;
        mat dataset;
        double fs = 0.0;
        load_dataset(data_dir + "/dataset_" + std::to_string(idx) + ".txt", dataset, fs);
        const int N = static_cast<int>(dataset[0].size());

        const double fD_res = fs / (2.0 * N);
        const double range_res = (3.0 * 100000000.0) / fs;
        const int Doppler_cells = static_cast<int>(max_Doppler / fD_res) + search_window_size + win[0] + 1;
        const int max_Doppler_ext = static_cast<int>(std::ceil(Doppler_cells * fD_res));
        const int max_range_ext = static_cast<int>(std::pow(2, std::ceil(std::log2((max_range / range_res) + search_window_size + win[0]))));

        const int ref_target_range_cell = static_cast<int>(std::round(target_reference_track[t][7] / range_res));
        const int ref_target_Doppler_cell = static_cast<int>(std::round(target_reference_track[t][8] / fD_res) + Doppler_cells);
        double ref_target_doa_dir = target_reference_track[t][9];
        if (ref_target_doa_dir < 0) ref_target_doa_dir += 180;
        if (ref_target_doa_dir > 180) ref_target_doa_dir -= 180;

        const std::string dir = out_root + "/dump/" + std::to_string(idx);
        std::filesystem::create_directories(dir);

        // Stage 1: Channel isolation
        std::pair<array, mat> ref_and_survs = isolate_channels(dataset, ref_ch_index);
        array& ref_ch = ref_and_survs.first;
        mat& surv_chs = ref_and_survs.second;
        dump_1d(dir + "/01_ref_ch.bin", ref_ch);
        for (size_t m = 0; m < surv_chs.size(); ++m)
            dump_1d(dir + "/02_surv_ch_" + std::to_string(m) + ".bin", surv_chs[m]);

        mat surv_chs_t = transpose(surv_chs);

        // Stage 2: Beamforming (replicate estimate_corr_matrix + optimal_Wiener_beamform + normalize)
        const std::complex<double> ci(0, 2 * M_PI);
        mat R = estimate_corr_matrix(surv_chs_t);
        mat aS(array_alignment.size(), array(1));
        for (size_t j = 0; j < array_alignment.size(); ++j)
            aS[j][0] = std::exp(array_alignment[j] * ci * std::cos(0.0174533 * ref_target_doa_dir));
        mat w = optimal_Wiener_beamform(R, aS);
        double sum_w = 0;
        for (size_t j = 0; j < w.size(); ++j)
            for (size_t k = 0; k < w[j].size(); ++k)
                sum_w += std::abs(w[j][k]);
        const double norm_factor = std::sqrt(static_cast<double>(M)) / sum_w;
        for (size_t j = 0; j < w.size(); ++j)
            for (size_t k = 0; k < w[j].size(); ++k)
                w[j][k] *= norm_factor;
        dump_weights(dir + "/03_beam_w.bin", w);

        array surv_ch = beamform_surveillance(surv_chs_t, array_alignment, ref_target_doa_dir);
        dump_1d(dir + "/03_beam_surv.bin", surv_ch);

        // Stage 3: Clutter cancellation
        array filt_surv_ch = time_domain_filter_surveillance(ref_ch, surv_ch, K);
        dump_1d(dir + "/04_clutter.bin", filt_surv_ch);

        // Stage 4: Windowing
        array windowed_surv_ch = windowing(filt_surv_ch, cc_det_window);
        dump_1d(dir + "/05_windowed.bin", windowed_surv_ch);

        // Stage 5: Detection - dump RAW rd_matrix BEFORE find_target abs-ifies it
        mat rd_matrix = cc_detector_ons(ref_ch, windowed_surv_ch, static_cast<int>(fs), max_Doppler_ext, max_range_ext);
        dump_2d(dir + "/06_rdmap.bin", rd_matrix);

        // Stage 6: Target finding
        std::pair<int, int> results = find_target(rd_matrix,
            ref_target_range_cell - search_window_size,
            ref_target_range_cell + search_window_size,
            ref_target_Doppler_cell - search_window_size,
            ref_target_Doppler_cell + search_window_size);
        const int target_range_cell = results.first;
        const int target_doppler_cell = results.second;

        // Stage 7: DOA per-channel maps + azimuth
        std::vector<std::pair<int, int>> hit_list = {{target_range_cell, target_doppler_cell}};
        tens rd_maps_for_doa(surv_chs.size());
        for (size_t m = 0; m < surv_chs.size(); ++m) {
            array f = time_domain_filter_surveillance(ref_ch, surv_chs[m], K);
            array wnd = windowing(f, cc_det_window);
            rd_maps_for_doa[m] = cc_detector_ons(ref_ch, wnd, static_cast<int>(fs), max_Doppler_ext, max_range_ext);
            dump_2d(dir + "/07_doa_rd_" + std::to_string(m) + ".bin", rd_maps_for_doa[m]);
        }
        std::vector<double> estimated_doas = target_DOA_estimation(rd_maps_for_doa, hit_list, array_alignment);
        const double azimuth_est = estimated_doas.empty() ? 0.0 : estimated_doas[0];

        // Stage 8: Metric extraction
        std::vector<int> target_rd = {target_range_cell, target_doppler_cell};
        target_rd[1] -= (static_cast<int>(rd_matrix.size()) - 1) / 2;
        win_pos[1] -= (static_cast<int>(rd_matrix.size()) - 1) / 2;
        const double delta = extract_delta(rd_matrix, target_rd, win);

        dump_results(dir + "/08_results.bin",
                     static_cast<double>(ref_target_range_cell), static_cast<double>(target_range_cell),
                     static_cast<double>(ref_target_Doppler_cell), static_cast<double>(target_doppler_cell),
                     target_reference_track[t][9], azimuth_est, delta);

        std::cout << "idx=" << idx << " range " << ref_target_range_cell << "/" << target_range_cell
                  << " doppler " << ref_target_Doppler_cell << "/" << target_doppler_cell
                  << " az " << target_reference_track[t][9] << "/" << azimuth_est
                  << " sinr " << delta << std::endl;
    }
    return 0;
}