#ifndef HAMPR_IO_LIVE_DATA_SOURCE_HPP
#define HAMPR_IO_LIVE_DATA_SOURCE_HPP

#include <hampr/io/data_source.hpp>
#include <hampr/io/buffer_pool.hpp>
#include <hampr/io/radio_source.hpp>
#include <memory>

namespace hampr {

// LiveDataSource: implements DataSource using a RadioSource.
// No SDR coupling in DSP — this is the only SDR-aware component.
class LiveDataSource : public DataSource {
public:
    explicit LiveDataSource(std::unique_ptr<RadioSource> source);

    // DataSource interface
    IQMatrix load(const std::string& filename, double& fs) override;
    TargetTrack load_track(const std::string& filename) override;

    // Streaming interface (via RadioSource + BufferPool)
    bool open(const std::string& device_spec);
    IQMatrix next_batch(size_t batch_size);
    bool has_more() const;
    void reset();
    void close();
    double sampling_rate() const;
    int num_channels() const;
    int total_samples() const;
    int current_position() const;

    // RadioSource configuration passthrough
    void set_frequency(double freq_hz) { source_->set_frequency(freq_hz); }
    void set_gain(double gain_db) { source_->set_gain(gain_db); }
    void set_sample_rate(double fs) { source_->set_sample_rate(fs); }

private:
    std::unique_ptr<RadioSource> source_;
    std::unique_ptr<BufferPool> buffer_pool_;
};

} // namespace hampr

#endif // HAMPR_IO_LIVE_DATA_SOURCE_HPP
