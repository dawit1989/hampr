#ifndef HAMPR_IO_CONSOLE_SINK_HPP
#define HAMPR_IO_CONSOLE_SINK_HPP

#include <hampr/io/data_sink.hpp>

namespace hampr {

class ConsoleSink : public DataSink {
public:
    void write(const ProcessingResult& result) override;
    void write_summary(const std::vector<ProcessingResult>& results) override;
};

} // namespace hampr

#endif // HAMPR_IO_CONSOLE_SINK_HPP
