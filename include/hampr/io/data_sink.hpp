#ifndef HAMPR_IO_DATA_SINK_HPP
#define HAMPR_IO_DATA_SINK_HPP

#include <hampr/core/types.hpp>
#include <string>

namespace hampr {

class DataSink {
public:
    virtual ~DataSink() = default;
    virtual void write(const ProcessingResult& result) = 0;
    virtual void write_summary(const std::vector<ProcessingResult>& results) = 0;
};

} // namespace hampr

#endif // HAMPR_IO_DATA_SINK_HPP
