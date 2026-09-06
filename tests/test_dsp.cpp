#include <hampr/dsp/channel_isolator.hpp>
#include <hampr/dsp/detector.hpp>
#include <hampr/dsp/metric_extractor.hpp>
#include <hampr/core/config.hpp>
#include <hampr/dsp/pipeline.hpp>
#include <cassert>
#include <iostream>

int main() {
    using namespace hampr;

    // Test channel isolation
    {
        IQMatrix data = {
            {complex(1,0), complex(2,0), complex(3,0)},
            {complex(4,0), complex(5,0), complex(6,0)},
            {complex(7,0), complex(8,0), complex(9,0)},
            {complex(10,0), complex(11,0), complex(12,0)}
        };
        ChannelIsolator isolator(0);
        auto result = isolator.isolate(data);
        assert(result.first.size() == 3);
        assert(result.second.size() == 3);
        assert(result.first[0] == complex(1,0));
    }

    // Test Hanni window
    {
        Detector detector("Hann");
        auto window = detector.hann(4);
        assert(window.size() == 4);
        // Hann(0) in [0, M-1] should be 0
        assert(std::abs(window[0]) < 1e-10);
    }

    // Test windowing
    {
        Detector detector("Hann");
        IQBuffer data = {complex(1,0), complex(2,0), complex(3,0), complex(4,0)};
        auto windowed = detector.apply_window(data);
        assert(windowed.size() == 4);
    }

    // Test pipeline creation
    {
        Config config;
        Pipeline pipeline(config);
        // Just verify it constructs without error
        assert(true);
    }

    std::cout << "All DSP tests passed!" << std::endl;
    return 0;
}

