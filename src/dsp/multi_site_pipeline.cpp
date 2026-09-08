#include <hampr/dsp/multi_site_pipeline.hpp>
#include <hampr/core/exception.hpp>
#include <hampr/dsp/array_geometry.hpp>
#include <chrono>

namespace hampr {

static const double NANOSECONDS_IN_SECOND = 1000000000.0;

MultiSitePipeline::MultiSitePipeline(const Config& config)
    : config_(config),
      synchronizer_(std::make_unique<Synchronizer>()),
      fusion_(std::make_unique<ResultFusion>()),
      localizer_(std::make_unique<TDoALocalizer>()),
      fusion_engine_(std::make_unique<MultiStaticFusion>()),
      multi_static_(std::make_unique<MultiStaticCorrelator>()) {

    if (!config_.multi_site_receivers.empty()) {
        for (const auto& info : config_.multi_site_receivers) {
            auto site_config = config;
            site_config.ref_channel_index = info.ref_channel_index;
            site_config.num_antennas = info.num_channels - 1;
            site_config.antenna_spacing = info.array_geometry.element_spacing;
            site_pipelines_.push_back(std::make_unique<Pipeline>(site_config));
        }
    } else {
        site_pipelines_.push_back(std::make_unique<Pipeline>(config));
    }
}

ProcessingResult MultiSitePipeline::process_site(
        int site_idx,
        const IQMatrix& iq_data,
        const TargetTrackPoint& ref_track,
        double fs, int total_samples) {

    if (site_idx >= static_cast<int>(site_pipelines_.size()) || site_idx < 0)
        throw HamprException("MultiSitePipeline: site index out of range");

    return site_pipelines_[site_idx]->process(iq_data, ref_track, fs, total_samples);
}

MultiSiteResult MultiSitePipeline::build_result(
        const std::vector<ProcessingResult>& site_results,
        const std::vector<ReceiverInfo>& receivers,
        double processing_time) {

    MultiSiteResult result;
    result.site_results = site_results;
    result.total_time = processing_time;

    if (site_results.empty())
        return result;

    if (config_.enable_fusion && site_results.size() >= 2) {
        auto fused = fusion_->fuse_bearings(site_results, receivers);
        result.latitude = fused.latitude;
        result.longitude = fused.longitude;
        result.altitude = fused.altitude;
    }

    return result;
}

MultiSiteResult MultiSitePipeline::process(const MultiSiteData& data,
                                            const TargetTrack& track) {
    auto t_start = std::chrono::high_resolution_clock::now();

    if (data.iq_data.empty() || track.empty()) {
        MultiSiteResult empty;
        empty.total_time = 0.0;
        return empty;
    }

    auto t_sync_start = std::chrono::high_resolution_clock::now();
    SynchronizedBatch batch = synchronizer_->synchronize(data);
    double sync_time = (std::chrono::high_resolution_clock::now() - t_sync_start).count() / NANOSECONDS_IN_SECOND;

    auto t_fuse_start = std::chrono::high_resolution_clock::now();
    MultiSiteResult result = process_batch(batch, track);
    double fuse_time = (std::chrono::high_resolution_clock::now() - t_fuse_start).count() / NANOSECONDS_IN_SECOND;

    result.fusion_time = sync_time + fuse_time;
    result.total_time = (std::chrono::high_resolution_clock::now() - t_start).count() / NANOSECONDS_IN_SECOND;
    return result;
}

MultiSiteResult MultiSitePipeline::process_batch(const SynchronizedBatch& batch,
                                                  const TargetTrack& track) {
    if (batch.aligned_data.empty() || track.empty())
        return MultiSiteResult{};

    const auto& ref_track = track.front();
    double fs = batch.fs;
    int num_samples = batch.batch_samples;

    std::vector<ProcessingResult> site_results;
    site_times_.clear();

    for (size_t s = 0; s < batch.aligned_data.size(); ++s) {
        auto t0 = std::chrono::high_resolution_clock::now();

        ProcessingResult site_result = process_site(
            static_cast<int>(s), batch.aligned_data[s], ref_track, fs, num_samples);

        site_times_.push_back(
            (std::chrono::high_resolution_clock::now() - t0).count() / NANOSECONDS_IN_SECOND);

        site_results.push_back(site_result);
    }

    TargetTrack single_track = {ref_track};

    // Multi-static cross-correlation
    auto ms_results = multi_static_->cross_correlate(batch);
    auto bistatic_ranges = multi_static_->compute_bistatic_ranges(ms_results, batch.fs);
    double ms_confidence = multi_static_->cooperative_confidence(ms_results, batch.aligned_data[0]);

    MultiSiteResult result = build_result(site_results, batch.receivers, 0.0);
    result.ms_confidence = ms_confidence;

    if (config_.enable_localization && batch.aligned_data.size() >= 2) {
        auto geo = localizer_->localize(batch, single_track);
        result.latitude = geo.latitude;
        result.longitude = geo.longitude;
        result.altitude = geo.altitude;
    }

    return result;
}

} // namespace hampr
