#include <hampr/io/streaming_source.hpp>
#include <iostream>

namespace hampr {

StreamingDataSource::StreamingDataSource(size_t buffer_capacity)
    : buffer_capacity_(buffer_capacity) {}

IQMatrix StreamingDataSource::load(const std::string& filename, double& fs) {
    open(filename);
    fs = fs_;
    int N = total_samples_;
    IQMatrix data(num_channels_, IQBuffer(N));
    for (int i = 0; i < num_channels_; ++i)
        for (int j = 0; j < N; ++j)
            data[i][j] = channels_[i][j];
    return data;
}

TargetTrack StreamingDataSource::load_track(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open())
        throw HamprException("Cannot open target reference track: " + filename);

    int rows, cols;
    infile >> rows >> cols;

    TargetTrack track;
    track.reserve(rows);
    for (int i = 0; i < rows; ++i) {
        TargetTrackPoint pt;
        double time_idx;
        infile >> time_idx >> pt.timestamp >> pt.latitude >> pt.longitude
               >> pt.altitude >> pt.speed >> pt.direction >> pt.range
               >> pt.doppler >> pt.azimuth;
        pt.time_index = static_cast<int>(time_idx);
        track.push_back(pt);
    }
    return track;
}

bool StreamingDataSource::open(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open())
        throw HamprException("Cannot open dataset file: " + filename);

    infile >> fs_;

    int rows, cols;
    infile >> rows >> cols;

    num_channels_ = rows;
    total_samples_ = cols;

    channels_.resize(rows, IQBuffer(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) {
            double real, imag;
            infile >> real >> imag;
            channels_[i][j] = complex(real, imag);
        }

    buffer_pool_ = std::make_unique<BufferPool>(buffer_capacity_);
    current_sample_ = 0;
    is_open_ = true;

    return true;
}

IQMatrix StreamingDataSource::next_batch(size_t batch_size) {
    if (!is_open_)
        throw HamprException("StreamingDataSource is not open");

    int remaining = total_samples_ - current_sample_;
    int actual_batch = static_cast<int>(std::min(batch_size, static_cast<size_t>(remaining)));

    IQMatrix batch(num_channels_, IQBuffer(actual_batch));
    for (int ch = 0; ch < num_channels_; ++ch)
        for (int i = 0; i < actual_batch; ++i)
            batch[ch][i] = channels_[ch][current_sample_ + i];

    current_sample_ += actual_batch;
    return batch;
}

bool StreamingDataSource::has_more() const {
    return is_open_ && current_sample_ < total_samples_;
}

void StreamingDataSource::reset() {
    current_sample_ = 0;
}

void StreamingDataSource::close() {
    is_open_ = false;
    buffer_pool_.reset();
}

} // namespace hampr
