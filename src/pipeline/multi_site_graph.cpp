#include <hampr/pipeline/multi_site_graph.hpp>
#include <iostream>

namespace hampr {

static const double NANOSECONDS_IN_SECOND = 1000000000.0;
static const double MILLISECONDS_IN_SECOND = 1000.0;

MultiSiteGraph::MultiSiteGraph(const Config& config)
    : pipeline_(std::make_unique<MultiSitePipeline>(config)),
      synchronizer_(std::make_unique<Synchronizer>()),
      start_time_(std::chrono::steady_clock::now()) {}

std::vector<MultiSiteResult> MultiSiteGraph::process_stream(
        MultiSiteDataSource& source,
        const TargetTrack& track,
        size_t batch_size) {

    std::vector<MultiSiteResult> results;
    if (track.empty())
        return results;

    const auto& ref_track = track.front();
    source.reset_all();
    size_t batch_idx = 0;

    while (true) {
        // Check if any site has more data
        bool has_any = false;
        std::vector<IQMatrix> site_batches;
        std::vector<ReceiverInfo> active_receivers;

        // Build per-site batches
        for (size_t s = 0; s < source.receivers().size(); ++s) {
            const auto& info = source.receivers()[s];
            if (!source.has_more_site(info.id))
                continue;
            has_any = true;
            IQMatrix batch = source.next_batch_site(info.id, batch_size);
            if (!batch.empty())
                site_batches.push_back(batch);
            active_receivers.push_back(info);
        }

        if (!has_any || site_batches.empty())
            break;

        auto batch_start = std::chrono::high_resolution_clock::now();

        // Build MultiSiteData from batches
        MultiSiteData data;
        data.fs = source.sampling_rate();
        data.receivers = active_receivers;
        data.iq_data = site_batches;

        // Synchronize
        SynchronizedBatch sync_batch = synchronizer_->synchronize(data);

        // Process
        MultiSiteResult result = pipeline_->process_batch(sync_batch, track);

        auto batch_end = std::chrono::high_resolution_clock::now();
        double latency_ms = (batch_end - batch_start).count() /
                            (NANOSECONDS_IN_SECOND / MILLISECONDS_IN_SECOND);

        {
            std::lock_guard<std::mutex> lock(results_mtx_);
            results.push_back(result);
            latencies_ms_.push_back(latency_ms);
            batches_processed_++;
            if (!sync_batch.aligned_data.empty() && !sync_batch.aligned_data[0].empty())
                total_samples_processed_ += sync_batch.batch_samples;
        }

        std::cout << "Batch " << batch_idx << ": " << sync_batch.batch_samples
                  << " samples, latency=" << result.total_time << "s"
                  << ", wall=" << latency_ms << "ms" << std::endl;
        batch_idx++;
    }

    return results;
}

std::future<MultiSiteResult> MultiSiteGraph::process_batch_async(
        const SynchronizedBatch& batch,
        const TargetTrack& track) {

    auto future = std::async(std::launch::async, [this, batch, track]() {
        auto t_start = std::chrono::high_resolution_clock::now();
        MultiSiteResult result = pipeline_->process_batch(batch, track);
        auto t_end = std::chrono::high_resolution_clock::now();
        double latency_ms = (t_end - t_start).count() /
                            (NANOSECONDS_IN_SECOND / MILLISECONDS_IN_SECOND);

        std::lock_guard<std::mutex> lock(results_mtx_);
        latencies_ms_.push_back(latency_ms);
        batches_processed_++;
        if (!batch.aligned_data.empty() && !batch.aligned_data[0].empty())
            total_samples_processed_ += batch.batch_samples;
        return result;
    });
    futures_.push_back(std::move(future));
    return std::move(futures_.back());
}

void MultiSiteGraph::wait_all() {
    for (auto& f : futures_)
        f.wait();
    futures_.clear();
}

double MultiSiteGraph::avg_latency_ms() const {
    if (latencies_ms_.empty())
        return 0.0;
    double sum = 0.0;
    for (double lat : latencies_ms_)
        sum += lat;
    return sum / latencies_ms_.size();
}

double MultiSiteGraph::throughput_samples_per_sec() const {
    auto elapsed = std::chrono::steady_clock::now() - start_time_;
    double elapsed_seconds = elapsed.count() / NANOSECONDS_IN_SECOND;
    if (elapsed_seconds <= 0)
        return 0.0;
    return total_samples_processed_ / elapsed_seconds;
}

} // namespace hampr
