#include <hampr/dsp/pipeline.hpp>
#include <hampr/io/text_data_source.hpp>
#include <hampr/io/console_sink.hpp>
#include <hampr/io/file_sink.hpp>
#include <hampr/core/config.hpp>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    std::string data_dir = "dataset";
    if (argc > 1) data_dir = argv[1];

    hampr::Config config;
    config.ref_channel_index = 0;
    config.filter_taps = 128;
    config.num_antennas = 3;
    config.antenna_spacing = 0.528;
    config.window_type = "Hann";
    config.search_window_size = 8;
    config.metric_win = {6, 6, 3, 3};

    hampr::Pipeline pipeline(config);
    hampr::TextDataSource data_source;
    hampr::ConsoleSink console_sink;
    hampr::FileSink file_sink("results.txt");

    auto track = data_source.load_track(data_dir + "/target_reference_track.txt");

    int start_idx = track.front().time_index;
    int stop_idx = track.back().time_index;

    std::vector<hampr::ProcessingResult> results;

    for (int t = 0; t <= stop_idx - start_idx; ++t) {
        std::cout << "====================== Processing Time Index "
                  << t + start_idx << " ======================\n\n";

        double fs;
        std::string filename = data_dir + "/dataset_" + std::to_string(t + start_idx) + ".txt";
        auto iq_data = data_source.load(filename, fs);
        int N = static_cast<int>(iq_data[0].size());

        auto result = pipeline.process(iq_data, track[t], fs, N);

        console_sink.write(result);
        file_sink.write(result);
        results.push_back(result);
    }

    console_sink.write_summary(results);

    return 0;
}
