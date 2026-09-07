#ifndef HAMPR_IO_STREAMING_SOURCE_HPP
#define HAMPR_IO_STREAMING_SOURCE_HPP

#include <hampr/io/data_source.hpp>
#include <hampr/io/buffer_pool.hpp>
#include <hampr/core/exception.hpp>
#include <fstream>
#include <memory>

namespace hampr {

// Streaming data source that provides IQ data in batches.
// Supports file replay as a deterministic test source.
class StreamingDataSource : public DataSource {
public:
    StreamingDataSource(size_t buffer_capacity = 65536);

    // File-based DataSource interface (loads entire file)
    IQMatrix load(const std::string& filename, double& fs) override;
    TargetTrack load_track(const std::string& filename) override;

    // Streaming interface
    bool open(const std::string& filename);
    IQMatrix next_batch(size_t batch_size);
    bool has_more() const;
    void reset();
    void close();
    double sampling_rate() const { return fs_; }
    int num_channels() const { return num_channels_; }
    int total_samples() const { return total_samples_; }
    int current_position() const { return current_sample_; }

private:
    size_t buffer_capacity_;
    std::unique_ptr<BufferPool> buffer_pool_;
    std::vector<IQBuffer> channels_;
    int current_sample_ = 0;
    int total_samples_ = 0;
    int num_channels_ = 0;
    double fs_ = 0.0;
    bool is_open_ = false;
};

} // namespace hampr

#endif // HAMPR_IO_STREAMING_SOURCE_HPP
