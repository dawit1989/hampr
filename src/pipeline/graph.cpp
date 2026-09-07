#include <hampr/pipeline/graph.hpp>
#include <iostream>

namespace hampr {

static const double NANOSECONDS_IN_SECOND = 1000000000.0;
static const double MILLISECONDS_IN_SECOND = 1000.0;

PipelineGraph::PipelineGraph(const Config& config)
    : pipeline_(std::make_unique<Pipeline>(config)),
      start_time_(std::chrono::steady_clock::now()) {}

ProcessingResult PipelineGraph::process_batch(
        const IQMatrix& iq_data,
        const TargetTrackPoint& ref_track,
        double fs,
        int total_samples) {

    auto t_start = std::chrono::high_resolution_clock::now();

    ProcessingResult result = pipeline_->process(iq_data, ref_track, fs, total_samples);

    auto t_end = std::chrono::high_resolution_clock::now();
    double latency_ms = (t_end - t_start).count() /
                        (NANOSECONDS_IN_SECOND / MILLISECONDS_IN_SECOND);

    {
        std::lock_guard<std::mutex> lock(results_mtx_);
        results_.push_back(result);
        latencies_ms_.push_back(latency_ms);
        batches_processed_++;
        total_samples_processed_ += total_samples;
    }

    return result;
}

std::vector<ProcessingResult> PipelineGraph::process_stream(
        StreamingDataSource& source,
        const TargetTrack& track,
        size_t batch_size) {

    std::vector<ProcessingResult> results;
    if (track.empty())
        return results;

    const auto& ref_track = track.front();
    double fs = source.sampling_rate();

    source.reset();
    size_t batch_idx = 0;

    while (source.has_more()) {
        auto batch_start = std::chrono::high_resolution_clock::now();

        IQMatrix batch = source.next_batch(batch_size);
        int batch_samples = static_cast<int>(batch[0].size());
        if (batch_samples == 0)
            break;

        TargetTrackPoint track_point = ref_track;
        track_point.range = static_cast<double>(batch_samples);

        auto result = process_batch(batch, track_point, fs, batch_samples);
        results.push_back(result);

        auto batch_end = std::chrono::high_resolution_clock::now();
        double elapsed_ms = (batch_end - batch_start).count() /
                            (NANOSECONDS_IN_SECOND / MILLISECONDS_IN_SECOND);

        std::cout << "Batch " << batch_idx << ": " << batch_samples
                  << " samples, latency=" << result.total_time << "s"
                  << ", wall=" << elapsed_ms << "ms" << std::endl;
        batch_idx++;
    }

    return results;
}

std::future<ProcessingResult> PipelineGraph::process_batch_async(
        const IQMatrix& iq_data,
        const TargetTrackPoint& ref_track,
        double fs,
        int total_samples) {

    auto future = std::async(std::launch::async, [this, iq_data, ref_track, fs, total_samples]() {
        return process_batch(iq_data, ref_track, fs, total_samples);
    });
    futures_.push_back(std::move(future));
    return std::move(futures_.back());
}

void PipelineGraph::wait_all() {
    for (auto& f : futures_) {
        f.wait();
    }
    futures_.clear();
}

double PipelineGraph::avg_latency_ms() const {
    if (latencies_ms_.empty())
        return 0.0;
    double sum = 0.0;
    for (double lat : latencies_ms_)
        sum += lat;
    return sum / latencies_ms_.size();
}

double PipelineGraph::throughput_samples_per_sec() const {
    auto elapsed = std::chrono::steady_clock::now() - start_time_;
    double elapsed_seconds = elapsed.count() / NANOSECONDS_IN_SECOND;
    if (elapsed_seconds <= 0)
        return 0.0;
    return total_samples_processed_ / elapsed_seconds;
}

} // namespace hampr
