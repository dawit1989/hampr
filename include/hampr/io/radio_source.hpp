#ifndef HAMPR_IO_RADIO_SOURCE_HPP
#define HAMPR_IO_RADIO_SOURCE_HPP

#include <hampr/core/types.hpp>
#include <string>
#include <vector>
#include <memory>

namespace hampr {

// Abstract interface for SDR hardware backends.
// Implementations: FileRadioSource (mock), USDRadioSource (UHD), RTLRadioSource (RTL-SDR)
class RadioSource {
public:
    virtual ~RadioSource() = default;

    virtual bool open(const std::string& device) = 0;
    virtual void close() = 0;
    virtual bool start() = 0;
    virtual void stop() = 0;
    virtual bool is_open() const = 0;

    // Read IQ samples into buffer (returns number of samples read)
    virtual size_t read(IQBuffer& buffer, size_t max_samples) = 0;

    // Configuration
    virtual void set_frequency(double freq_hz) = 0;
    virtual void set_gain(double gain_db) = 0;
    virtual void set_sample_rate(double fs) = 0;
    virtual void set_center_frequency(double freq_hz) { set_frequency(freq_hz); }
    virtual void set_bandwidth(double bw_hz) {}

    virtual double sampling_rate() const = 0;
    virtual double frequency() const = 0;
    virtual double gain() const = 0;
    virtual int num_channels() const = 0;
};

// File-based mock SDR source for testing.
// Replays IQ data from text files without requiring hardware.
class FileRadioSource : public RadioSource {
public:
    FileRadioSource();
    ~FileRadioSource() override = default;

    bool open(const std::string& filename) override;
    void close() override;
    bool start() override;
    void stop() override;
    bool is_open() const override;

    size_t read(IQBuffer& buffer, size_t max_samples) override;

    void set_frequency(double freq_hz) override { freq_hz_ = freq_hz; }
    void set_gain(double gain_db) override { gain_db_ = gain_db; }
    void set_sample_rate(double fs) override { fs_ = fs; }
    double sampling_rate() const override { return fs_; }
    double frequency() const override { return freq_hz_; }
    double gain() const override { return gain_db_; }
    int num_channels() const override { return 1; }

    bool has_more() const;
    void reset();
    int total_samples() const { return static_cast<int>(iq_data_.size()); }
    int current_position() const { return position_; }

private:
    IQBuffer iq_data_;
    int position_ = 0;
    double fs_ = 0.0;
    double freq_hz_ = 0.0;
    double gain_db_ = 0.0;
    bool is_open_ = false;
};

} // namespace hampr

#endif // HAMPR_IO_RADIO_SOURCE_HPP
