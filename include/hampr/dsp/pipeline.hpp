#ifndef HAMPR_DSP_PIPELINE_HPP
#define HAMPR_DSP_PIPELINE_HPP

#include <hampr/core/types.hpp>
#include <hampr/core/config.hpp>
#include <hampr/dsp/channel_isolator.hpp>
#include <hampr/dsp/beamformer.hpp>
#include <hampr/dsp/clutter_canceler.hpp>
#include <hampr/dsp/detector.hpp>
#include <hampr/dsp/doa_estimator.hpp>
#include <hampr/dsp/metric_extractor.hpp>
#include <memory>

namespace hampr {

class Pipeline {
public:
    explicit Pipeline(const Config& config = Config());

    ProcessingResult process(const IQMatrix& iq_data,
                             const TargetTrackPoint& ref_track,
                             double fs,
                             int total_samples);

private:
    Config config_;
    std::unique_ptr<ChannelIsolator> isolator_;
    std::unique_ptr<Beamformer> beamformer_;
    std::unique_ptr<ClutterCanceler> canceler_;
    std::unique_ptr<Detector> detector_;
    std::unique_ptr<DOAEstimator> doa_estimator_;
    std::unique_ptr<MetricExtractor> metric_extractor_;
};

} // namespace hampr

#endif // HAMPR_DSP_PIPELINE_HPP
