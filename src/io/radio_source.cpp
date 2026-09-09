#include <hampr/io/radio_source.hpp>
#include <hampr/io/text_data_source.hpp>
#include <fstream>
#include <sstream>

namespace hampr {

FileRadioSource::FileRadioSource() = default;

bool FileRadioSource::open(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open())
        return false;

    iq_data_.clear();
    position_ = 0;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        double re, im;
        char comma;
        if (ss >> re >> comma >> im) {
            iq_data_.push_back(complex(re, im));
        } else {
            // Try single value (real only)
            ss.clear();
            ss.str(line);
            double val;
            if (ss >> val)
                iq_data_.push_back(complex(val, 0.0));
        }
    }

    if (iq_data_.empty())
        return false;

    is_open_ = true;
    return true;
}

void FileRadioSource::close() {
    iq_data_.clear();
    position_ = 0;
    is_open_ = false;
}

bool FileRadioSource::start() {
    return is_open_;
}

void FileRadioSource::stop() {
    // No-op for file source
}

bool FileRadioSource::is_open() const {
    return is_open_;
}

size_t FileRadioSource::read(IQBuffer& buffer, size_t max_samples) {
    if (!is_open_ || position_ >= static_cast<int>(iq_data_.size()))
        return 0;

    size_t count = std::min(max_samples, static_cast<size_t>(iq_data_.size() - position_));
    buffer.resize(count);
    for (size_t i = 0; i < count; ++i)
        buffer[i] = iq_data_[position_ + i];
    position_ += static_cast<int>(count);
    return count;
}

bool FileRadioSource::has_more() const {
    return is_open_ && position_ < static_cast<int>(iq_data_.size());
}

void FileRadioSource::reset() {
    position_ = 0;
}

} // namespace hampr
