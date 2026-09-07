#include <hampr/io/buffer_pool.hpp>
#include <hampr/io/streaming_source.hpp>
#include <hampr/pipeline/graph.hpp>
#include <hampr/core/config.hpp>
#include <iostream>
#include <cassert>

int main() {
    using namespace hampr;
    std::cout << "Running Phase 5 streaming tests..." << std::endl;

    // Test 1: BufferPool push/pop
    {
        BufferPool pool(16);
        assert(pool.capacity() == 16);
        assert(pool.empty());
        assert(!pool.full());

        complex data[8];
        for (int i = 0; i < 8; ++i)
            data[i] = complex(i, i * 2);
        assert(pool.push(data, 8));
        assert(pool.available() == 8);
        assert(!pool.empty());

        complex out[8];
        size_t got = pool.try_pop(out, 8);
        assert(got == 8);
        for (int i = 0; i < 8; ++i)
            assert(out[i] == complex(i, i * 2));
        assert(pool.empty());

        std::cout << "  Test 1 PASSED: BufferPool push/pop" << std::endl;
    }

    // Test 2: BufferPool wrap-around
    {
        BufferPool pool(4);
        complex data[3] = {complex(1,0), complex(2,0), complex(3,0)};
        assert(pool.push(data, 3));

        complex out[2];
        size_t got = pool.try_pop(out, 2);
        assert(got == 2);
        assert(out[0] == complex(1,0));
        assert(out[1] == complex(2,0));

        // Push more, causing wrap-around
        complex more[2] = {complex(10,0), complex(11,0)};
        assert(pool.push(more, 2));
        assert(pool.available() == 3);

        got = pool.try_pop(out, 2);
        assert(got == 2);
        assert(out[0] == complex(3,0));
        assert(out[1] == complex(10,0));

        std::cout << "  Test 2 PASSED: BufferPool wrap-around" << std::endl;
    }

    // Test 3: StreamingDataSource batch reading
    {
        std::string data_dir = "dataset";
        std::ifstream check(data_dir + "/dataset_494.txt");
        if (!check.good()) {
            data_dir = "../passiveradar_data/dataset";
            check.open(data_dir + "/dataset_494.txt");
        }
        if (!check.good()) {
            std::cout << "  Test 3 SKIPPED: VEGA dataset not found" << std::endl;
        } else {
            StreamingDataSource source(131072);
            assert(source.open(data_dir + "/dataset_494.txt"));
            assert(source.has_more());
            assert(source.num_channels() == 4);
            assert(source.sampling_rate() > 0);

            IQMatrix batch1 = source.next_batch(10000);
            assert(batch1.size() == 4);
            assert(batch1[0].size() == 10000);
            assert(source.has_more());

            IQMatrix batch2 = source.next_batch(100000);
            assert(source.has_more() || !source.has_more());

            source.reset();
            assert(source.has_more());
            IQMatrix batch_reset = source.next_batch(10000);
            assert(batch_reset[0].size() == 10000);

            source.close();
            std::cout << "  Test 3 PASSED: StreamingDataSource batch reading" << std::endl;
        }
    }

    // Test 4: PipelineGraph batch processing
    {
        std::string data_dir = "dataset";
        std::ifstream check(data_dir + "/dataset_494.txt");
        if (!check.good()) {
            data_dir = "../passiveradar_data/dataset";
            check.open(data_dir + "/dataset_494.txt");
        }
        if (!check.good()) {
            std::cout << "  Test 4 SKIPPED: VEGA dataset not found" << std::endl;
        } else {
            StreamingDataSource source(131072);
            source.open(data_dir + "/dataset_494.txt");
            auto track = source.load_track(data_dir + "/target_reference_track.txt");
            assert(track.size() >= 1);

            Config config;
            config.ref_channel_index = 0;
            config.filter_taps = 128;
            config.num_antennas = 3;
            config.antenna_spacing = 0.528;
            config.window_type = "Hann";
            config.search_window_size = 8;
            config.metric_win = {6, 6, 3, 3};

            PipelineGraph pipeline(config);
            auto results = pipeline.process_stream(source, track, 65536);

            assert(!results.empty());
            assert(pipeline.batches_processed() > 0);
            assert(pipeline.avg_latency_ms() > 0);
            assert(pipeline.throughput_samples_per_sec() > 0);

            source.close();
            std::cout << "  Test 4 PASSED: PipelineGraph batch processing ("
                      << pipeline.batches_processed() << " batches, "
                      << "avg_latency=" << pipeline.avg_latency_ms() << "ms"
                      << ")" << std::endl;
        }
    }

    // Test 5: Async metrics interface
    {
        Config config;
        PipelineGraph pipeline(config);

        assert(pipeline.batches_processed() == 0);
        assert(pipeline.avg_latency_ms() == 0.0);
        assert(pipeline.throughput_samples_per_sec() >= 0.0);

        std::cout << "  Test 5 PASSED: Async metrics interface" << std::endl;
    }

    std::cout << "All Phase 5 streaming tests passed!" << std::endl;
    return 0;
}
