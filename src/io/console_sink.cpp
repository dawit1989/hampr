#include <hampr/io/console_sink.hpp>
#include <iostream>
#include <iomanip>

namespace hampr {

void ConsoleSink::write(const ProcessingResult& result) {
    std::cout << "Target range cell ref/found ["
              << result.target_range_cell_ref << "/"
              << result.target_range_cell_found << "]\n";
    std::cout << "Target Doppler cell ref/found ["
              << result.target_doppler_cell_ref << "/"
              << result.target_doppler_cell_found << "]\n";
    std::cout << "Target azimuth ref/found ["
              << result.azimuth_ref << "/"
              << result.azimuth_estimated << "]\n";
    std::cout << "Target SINR " << result.sinr << "\n";

    std::cout << "\nTime to isolate channels: " << result.isolation_time << " seconds\n";
    std::cout << "Time for beamforming: " << result.beamform_time << " seconds\n";
    std::cout << "Time to perform time-domain filter surveillance: " << result.filter_time << " seconds\n";
    std::cout << "Time to detect: " << result.detector_time << " seconds\n";
    std::cout << "Time to find target: " << result.target_find_time << " seconds\n";
    std::cout << "Time to estimate direction of arrival: " << result.doa_time << " seconds\n";
    std::cout << "Total time to process: " << result.total_time << " seconds\n\n";
}

void ConsoleSink::write_summary(const std::vector<ProcessingResult>& results) {
    double avg_isolation = 0, avg_beamform = 0, avg_filter = 0;
    double avg_detect = 0, avg_find = 0, avg_doa = 0, avg_total = 0;

    for (const auto& r : results) {
        avg_isolation += r.isolation_time;
        avg_beamform += r.beamform_time;
        avg_filter += r.filter_time;
        avg_detect += r.detector_time;
        avg_find += r.target_find_time;
        avg_doa += r.doa_time;
        avg_total += r.total_time;
    }

    int n = static_cast<int>(results.size());
    std::cout << "====================== RESULTS SUMMARY ======================\n\n";
    std::cout << "Avg time to isolate channels: " << avg_isolation/n << " seconds\n";
    std::cout << "Avg time for beamforming: " << avg_beamform/n << " seconds\n";
    std::cout << "Avg time to perform time-domain filter surveillance: " << avg_filter/n << " seconds\n";
    std::cout << "Avg time to detect: " << avg_detect/n << " seconds\n";
    std::cout << "Avg time to find target: " << avg_find/n << " seconds\n";
    std::cout << "Avg time to estimate direction of arrival: " << avg_doa/n << " seconds\n";
    std::cout << "Avg total time to process single time index: " << avg_total/n << " seconds\n\n";
}

} // namespace hampr
