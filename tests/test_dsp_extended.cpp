#include <hampr/dsp/pipeline.hpp>
#include <hampr/dsp/channel_isolator.hpp>
#include <hampr/dsp/beamformer.hpp>
#include <hampr/dsp/clutter_canceler.hpp>
#include <hampr/dsp/detector.hpp>
#include <hampr/dsp/doa_estimator.hpp>
#include <hampr/core/exception.hpp>
#include <hampr/utils/fftw_wrapper.hpp>
#include <hampr/utils/math_utils.hpp>
#include <hampr/utils/linalg_utils.hpp>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <cassert>
#include <limits>

static double max_rel_error(const hampr::array& a, const hampr::array& b) {
    assert(a.size() == b.size());
    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double denom = std::abs(b[i]);
        double err = denom > 1e-15 ? std::abs(a[i] - b[i]) / denom : std::abs(a[i] - b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

int main() {
    using namespace hampr;
    std::cout << "Running extended DSP tests..." << std::endl;

    // Test 1: FFT/IFFT round-trip
    {
        array original(16);
        for (int i = 0; i < 16; ++i) original[i] = complex(i, i*0.5);
        array recovered = original;
        fft(recovered); ifft(recovered);
        assert(max_rel_error(original, recovered) < 1e-10);
        std::cout << "  Test 1 PASSED: FFT/IFFT round-trip" << std::endl;
    }

    // Test 2: ChannelIsolator throws on invalid ref_index
    {
        IQMatrix data(3, IQBuffer(3));
        ChannelIsolator isolator(5);
        bool threw = false;
        try { isolator.isolate(data); } catch (const HamprException&) { threw = true; }
        assert(threw);
        std::cout << "  Test 2 PASSED" << std::endl;
    }

    // Test 3: ChannelIsolator valid
    {
        IQMatrix data(4, IQBuffer(3));
        for (int i = 0; i < 4; ++i) data[i][0] = complex(i+1, 0);
        ChannelIsolator isolator(1);
        auto r = isolator.isolate(data);
        assert(r.first[0] == complex(2,0));
        assert(r.second[0][0] == complex(1,0));
        std::cout << "  Test 3 PASSED" << std::endl;
    }

    // Test 4: Beamformer mismatched dims throws
    {
        IQMatrix ant(2, IQBuffer(3));
        std::vector<double> alignment = {0.0, 0.528};
        Beamformer bf(0.528);
        bool threw = false;
        try { bf.beamform(ant, alignment, 90.0); } catch (const HamprException&) { threw = true; }
        assert(threw);
        std::cout << "  Test 4 PASSED" << std::endl;
    }

    // Test 5: ClutterCanceler mismatched sizes throws
    {
        ClutterCanceler canceler(128);
        IQBuffer ref(3), surv(2);
        bool threw = false;
        try { canceler.filter(ref, surv); } catch (const HamprException&) { threw = true; }
        assert(threw);
        std::cout << "  Test 5 PASSED" << std::endl;
    }

    // Test 6: DOAEstimator produces valid angle
    {
        DOAEstimator estimator;
        RDMap rd_map(5, IQBuffer(5, complex(0,0)));
        rd_map[2][2] = complex(1,0);
        std::vector<RDMap> rd_maps = {rd_map, rd_map, rd_map};
        std::vector<std::pair<int,int>> hit_list = {{2, 2}};
        std::vector<double> alignment = {0.0, 0.528, 1.056};
        auto doas = estimator.estimate(rd_maps, hit_list, alignment);
        assert(doas.size() == 1);
        assert(doas[0] >= 0.0 && doas[0] < 180.0);
        std::cout << "  Test 6 PASSED" << std::endl;
    }

    // Test 7: Hann window
    {
        Detector detector("Hann");
        auto w = detector.hann(8);
        assert(w.size() == 8);
        assert(std::abs(w[0]) < 1e-10);
        assert(std::abs(w[7]) < 1e-10);
        assert(std::abs(w[4] - 1.0) < 1e-10);
        std::cout << "  Test 7 PASSED" << std::endl;
    }

    // Test 8: MetricExtractor - matrix must be large enough
    {
        MetricExtractor extractor;
        int rows = 13;
        RDMap rd_matrix(rows, IQBuffer(rows, complex(0,0)));
        rd_matrix[10][4] = complex(10,0);
        std::vector<int> target_rd = {4, 4};
        std::vector<int> win = {2, 2, 1, 1};
        double snr = extractor.extract_snr(rd_matrix, target_rd, win);
        assert(std::isfinite(snr));
        std::cout << "  Test 8 PASSED: SNR=" << snr << std::endl;
    }

    std::cout << "All extended DSP tests passed!" << std::endl;
    return 0;
}
