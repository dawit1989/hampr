#ifndef HAMPR_DSP_MULTI_SITE_PIPELINE_HPP
#define HAMPR_DSP_MULTI_SITE_PIPELINE_HPP

#include <hampr/dsp/pipeline.hpp>
#include <hampr/dsp/synchronizer.hpp>
#include <hampr/dsp/fusion.hpp>
#include <hampr/dsp/localization.hpp>
#include <hampr/dsp/multi_static_correlator.hpp>
#include <hampr/dsp/multi_static_fusion.hpp>
#include <hampr/core/config.hpp>
#include <hampr/core/multi_site_types.hpp>
#include <memory>
#include <vector>

namespace hampr {

class MultiSitePipeline {
public:
    explicit MultiSitePipeline(const Config& config);

    MultiSiteResult process(const MultiSiteData& data,
                            const TargetTrack& track);

    MultiSiteResult process_batch(const SynchronizedBatch& batch,
                                  const TargetTrack& track);

    std::vector<double> site_processing_times() const { return site_times_; }

private:
    Config config_;
    std::unique_ptr<Synchronizer> synchronizer_;
    std::unique_ptr<ResultFusion> fusion_;
    std::unique_ptr<TDoALocalizer> localizer_;
    std::unique_ptr<MultiStaticFusion> fusion_engine_;
    std::unique_ptr<MultiStaticCorrelator> multi_static_;
    std::vector<std::unique_ptr<Pipeline>> site_pipelines_;
    std::vector<double> site_times_;

    ProcessingResult process_site(int site_idx,
                                  const IQMatrix& iq_data,
                                  const TargetTrackPoint& ref_track,
                                  double fs, int total_samples);

    MultiSiteResult build_result(const std::vector<ProcessingResult>& site_results,
                                 const std::vector<ReceiverInfo>& receivers,
                                 double processing_time);
};

} // namespace hampr

#endif // HAMPR_DSP_MULTI_SITE_PIPELINE_HPP
