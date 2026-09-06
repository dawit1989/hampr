#include <hampr/dsp/pipeline.hpp>
#include <hampr/utils/math_utils.hpp>
#include <chrono>
#include <cmath>

namespace hampr {

static const double NANOSECONDS_IN_SECOND = 1000000000.0;

static std::pair<int, int> find_target(RDMap& rd_matrix,
                                       int range_min, int range_max,
                                       int doppler_min, int doppler_max) {
    for (auto& row : rd_matrix)
        for (auto& val : row)
            val = complex(std::abs(val), 0);
    double max_val = 0;
    int argmax_d = 0, argmax_r = 0;
    for (int i = doppler_min; i < doppler_max; ++i)
        for (int j = range_min; j < range_max; ++j)
            if (rd_matrix[i][j].real() > max_val) {
                max_val = rd_matrix[i][j].real();
                argmax_d = i;
                argmax_r = j;
            }
    return {argmax_r, argmax_d};
}

Pipeline::Pipeline(const Config& config) : config_(config) {
    isolator_ = std::make_unique<ChannelIsolator>(config.ref_channel_index);
    beamformer_ = std::make_unique<Beamformer>(config.antenna_spacing);
    canceler_ = std::make_unique<ClutterCanceler>(config.filter_taps);
    detector_ = std::make_unique<Detector>(config.window_type);
    doa_estimator_ = std::make_unique<DOAEstimator>();
    metric_extractor_ = std::make_unique<MetricExtractor>();
}

ProcessingResult Pipeline::process(const IQMatrix& iq_data,
                                   const TargetTrackPoint& ref_track,
                                   double fs, int total_samples) {
    using namespace std::chrono;
    double max_range = ref_track.range;
    double max_doppler = ref_track.doppler;
    int N = total_samples;
    std::vector<double> array_alignment(config_.num_antennas);
    for (int i = 0; i < config_.num_antennas; ++i)
        array_alignment[i] = i * config_.antenna_spacing;
    double fD_res = fs / (2.0 * N);
    double range_res = 3.0e8 / fs;
    int doppler_cells = static_cast<int>(max_doppler / fD_res) + config_.search_window_size + config_.metric_win[0] + 1;
    int max_doppler_ext = static_cast<int>(std::ceil(doppler_cells * fD_res));
    int max_range_ext = static_cast<int>(std::pow(2, std::ceil(std::log2(max_range / range_res + config_.search_window_size + config_.metric_win[0]))));
    int ref_range_cell = static_cast<int>(std::round(ref_track.range / range_res));
    int ref_doppler_cell = static_cast<int>(std::round(ref_track.doppler / fD_res) + doppler_cells);
    double doa_dir = ref_track.azimuth;
    if (doa_dir < 0) doa_dir += 180;
    if (doa_dir > 180) doa_dir -= 180;
    ProcessingResult result;
    result.azimuth_ref = ref_track.azimuth;
    result.target_range_cell_ref = ref_range_cell;
    result.target_doppler_cell_ref = ref_doppler_cell;
    auto iter_start = high_resolution_clock::now();
    auto t0 = high_resolution_clock::now();
    auto ref_and_survs = isolator_->isolate(iq_data);
    auto& ref_ch = ref_and_survs.first;
    auto& surv_chs = ref_and_survs.second;
    auto surv_chs_t = transpose(surv_chs);
    result.isolation_time = (high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND;
    t0 = high_resolution_clock::now();
    auto beam_result = beamformer_->beamform(surv_chs_t, array_alignment, doa_dir);
    auto surv_ch = beam_result.surveillance_channel;
    result.beamform_time = (high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND;
    t0 = high_resolution_clock::now();
    auto filt_surv_ch = canceler_->filter(ref_ch, surv_ch);
    result.filter_time = (high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND;
    auto windowed_surv_ch = detector_->apply_window(filt_surv_ch);
    t0 = high_resolution_clock::now();
    RDMap rd_matrix = detector_->detect(ref_ch, windowed_surv_ch, fs, max_doppler_ext, max_range_ext);
    result.detector_time = (high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND;
    t0 = high_resolution_clock::now();
    auto target_result = find_target(rd_matrix,
        ref_range_cell - config_.search_window_size,
        ref_range_cell + config_.search_window_size,
        ref_doppler_cell - config_.search_window_size,
        ref_doppler_cell + config_.search_window_size);
    result.target_range_cell_found = target_result.first;
    result.target_doppler_cell_found = target_result.second;
    result.target_find_time = (high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND;
    t0 = high_resolution_clock::now();
    std::vector<std::pair<int,int>> hit_list = {{target_result.first, target_result.second}};
    std::vector<RDMap> rd_maps_for_doa(surv_chs.size());
    for (size_t m = 0; m < surv_chs.size(); ++m) {
        auto filt = canceler_->filter(ref_ch, surv_chs[m]);
        auto windowed = detector_->apply_window(filt);
        rd_maps_for_doa[m] = detector_->detect(ref_ch, windowed, fs, max_doppler_ext, max_range_ext);
    }
    auto doas = doa_estimator_->estimate(rd_maps_for_doa, hit_list, array_alignment);
    result.azimuth_estimated = doas[0];
    result.doa_time = (high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND;
    std::vector<int> win = config_.metric_win;
    std::vector<int> target_rd = {target_result.first, target_result.second};
    target_rd[1] -= (static_cast<int>(rd_matrix.size()) - 1) / 2;
    result.sinr = metric_extractor_->extract_snr(rd_matrix, target_rd, win);
    result.total_time = (high_resolution_clock::now() - iter_start).count() / NANOSECONDS_IN_SECOND;
    return result;
}

} // namespace hampr
