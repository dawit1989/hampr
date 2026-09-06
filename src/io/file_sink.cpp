#include <hampr/io/file_sink.hpp>
#include <fstream>
#include <hampr/core/exception.hpp>#include <iostream>

namespace hampr {

FileSink::FileSink(const std::string& filename)
    : filename_(filename) {}

FileSink::~FileSink() = default;

void FileSink::write(const ProcessingResult& result) {
    std::ofstream outfile(filename_, std::ios::app);
    if (!outfile.is_open()) {
        throw HamprException("Cannot open output file: " + filename_);
    }

    outfile << result.target_range_cell_ref << " "
            << result.target_range_cell_found << " "
            << result.target_doppler_cell_ref << " "
            << result.target_doppler_cell_found << " "
            << result.azimuth_ref << " "
            << result.azimuth_estimated << " "
            << result.sinr << "\n";
}

void FileSink::write_summary(const std::vector<ProcessingResult>& /*results*/) {
    // Summary is included in the per-result output
}

} // namespace hampr


