#ifndef HAMPR_IO_FILE_SINK_HPP
#define HAMPR_IO_FILE_SINK_HPP

#include <hampr/io/data_sink.hpp>
#include <string>

namespace hampr {

class FileSink : public DataSink {
public:
    explicit FileSink(const std::string& filename);
    ~FileSink();

    void write(const ProcessingResult& result) override;
    void write_summary(const std::vector<ProcessingResult>& results) override;

private:
    std::string filename_;
};

} // namespace hampr

#endif // HAMPR_IO_FILE_SINK_HPP
