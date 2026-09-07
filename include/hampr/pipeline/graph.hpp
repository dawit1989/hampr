#ifndef HAMPR_PIPELINE_GRAPH_HPP
#define HAMPR_PIPELINE_GRAPH_HPP

#include <hampr/dsp/pipeline.hpp>
#include <hampr/io/streaming_source.hpp>
#include <hampr/io/data_sink.hpp>
#include <hampr/core/config.hpp>
#include <vector>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <future>

namespace hampr {

// Pipeline graph scheduler for streaming/asynchronous batch processing.
// Measures per-batch latency and overall throughput.
class PipelineGraph {
public:
    explicit PipelineGraph(const Config& config = Config());

    // Process a single batch of data
    ProcessingResult process_batch(const IQMatrix& iq_data,
                                   const TargetTrackPoint& ref_track,
                                   double fs,
                                   int total_samples);

    // Process all data from a streaming source in batches
    std::vector<ProcessingResult> process_stream(StreamingDataSource& source,
                                                 const TargetTrack& track,
                                                 size_t batch_size);

    // Async processing: launch a batch for background processing
    std::future<ProcessingResult> process_batch_async(const IQMatrix& iq_data,
                                                      const TargetTrackPoint& ref_track,
                                                      double fs,
                                                      int total_samples);

    // Wait for all async batches to complete
    void wait_all();

    // Metrics
    double avg_latency_ms() const;
    double throughput_samples_per_sec() const;
    size_t batches_processed() const { return batches_processed_; }

private:
    std::unique_ptr<Pipeline> pipeline_;
    std::vector<ProcessingResult> results_;
    std::vector<double> latencies_ms_;
    std::vector<std::future<ProcessingResult>> futures_;

    std::atomic<size_t> batches_processed_{0};
    std::atomic<int> total_samples_processed_{0};
    std::chrono::steady_clock::time_point start_time_;

    mutable std::mutex results_mtx_;
};

} // namespace hampr

#endif // HAMPR_PIPELINE_GRAPH_HPP
