#include <hampr/io/live_data_source.hpp>
#include <hampr/io/text_data_source.hpp>

namespace hampr {

LiveDataSource::LiveDataSource(std::unique_ptr<RadioSource> source)
    : source_(std::move(source)),
      buffer_pool_(std::make_unique<BufferPool>(65536)) {}

bool LiveDataSource::open(const std::string& device_spec) {
    return source_->open(device_spec);
}

IQMatrix LiveDataSource::next_batch(size_t batch_size) {
    IQMatrix batch;
    batch.reserve(source_->num_channels());
    IQBuffer channel(batch_size);

    for (int ch = 0; ch < source_->num_channels(); ++ch) {
        size_t n = source_->read(channel, batch_size);
        channel.resize(n);
        batch.push_back(std::move(channel));
        channel = IQBuffer(batch_size);
    }

    return batch;
}

bool LiveDataSource::has_more() const {
    auto* file_src = dynamic_cast<FileRadioSource*>(source_.get());
    if (file_src)
        return file_src->has_more();
    return source_->is_open();
}

void LiveDataSource::reset() {
    // FileRadioSource has reset()
    auto* file_src = dynamic_cast<FileRadioSource*>(source_.get());
    if (file_src)
        file_src->reset();
}

void LiveDataSource::close() {
    source_->stop();
    source_->close();
}

double LiveDataSource::sampling_rate() const {
    return source_->sampling_rate();
}

int LiveDataSource::num_channels() const {
    return source_->num_channels();
}

int LiveDataSource::total_samples() const {
    auto* file_src = dynamic_cast<FileRadioSource*>(source_.get());
    return file_src ? file_src->total_samples() : 0;
}

int LiveDataSource::current_position() const {
    auto* file_src = dynamic_cast<FileRadioSource*>(source_.get());
    return file_src ? file_src->current_position() : 0;
}

IQMatrix LiveDataSource::load(const std::string& filename, double& fs) {
    TextDataSource loader;
    return loader.load(filename, fs);
}

TargetTrack LiveDataSource::load_track(const std::string& filename) {
    TextDataSource loader;
    return loader.load_track(filename);
}

} // namespace hampr
