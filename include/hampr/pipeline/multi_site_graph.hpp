#ifndef HAMPR_PIPELINE_MULTI_SITE_GRAPH_HPP
#define HAMPR_PIPELINE_MULTI_SITE_GRAPH_HPP

#include <hampr/dsp/multi_site_pipeline.hpp>
#include <hampr/io/multi_site_data_source.hpp>
#include <hampr/core/config.hpp>
#include <vector>
#include <atomic>
#include <chrono>
#include <mutex>
#include <future>
#include <thread>

namespace hampr {

class MultiSiteGraph {
public:
    explicit MultiSiteGraph(const Config& config);

    std::vector<MultiSiteResult> process_stream(
        MultiSiteDataSource& source,
        const TargetTrack& track,
        size_t batch_size);

    std::future<MultiSiteResult> process_batch_async(
        const SynchronizedBatch& batch,
        const TargetTrack& track);

    void wait_all();

    double avg_latency_ms() const;
    double throughput_samples_per_sec() const;
    size_t batches_processed() const { return batches_processed_; }

private:
    std::unique_ptr<MultiSitePipeline> pipeline_;
    std::unique_ptr<Synchronizer> synchronizer_;
    std::vector<MultiSiteResult> results_;
    std::vector<double> latencies_ms_;
    std::vector<std::future<MultiSiteResult>> futures_;

    std::atomic<size_t> batches_processed_{0};
    std::atomic<int> total_samples_processed_{0};
    std::chrono::steady_clock::time_point start_time_;

    mutable std::mutex results_mtx_;
};

} // namespace hampr

#endif // HAMPR_PIPELINE_MULTI_SITE_GRAPH_HPP
