#include <hampr/dsp/detector.hpp>
#include <cmath>
#include <cassert>

namespace hampr {

std::vector<double> Detector::hann(int M) const {
    std::vector<double> w(M);
    for (int n = 0; n < M; ++n)
        w[n] = 0.5 - 0.5 * std::cos(2.0 * M_PI * n / (M - 1));
    return w;
}

IQBuffer Detector::apply_window(const IQBuffer& surv_ch) const {
    std::vector<double> window;
    if (window_type_ == "Hann")
        window = hann(static_cast<int>(surv_ch.size()));

    int M = static_cast<int>(surv_ch.size());
    IQBuffer output(M);
    for (int i = 0; i < M; ++i)
        output[i] = surv_ch[i] * window[i];
    return output;
}

RDMap Detector::detect(IQBuffer& ref_ch, IQBuffer& surv_ch,
                       double fs, int fD_max, int r_max) const {
    return cc_detector_ons(ref_ch, surv_ch, fs, fD_max, r_max);
}

RDMap Detector::cc_detector_ons(IQBuffer& ref_ch, IQBuffer& surv_ch,
                                double fs, int fD_max, int r_max) const {
    int N = static_cast<int>(ref_ch.size());
    assert(N % r_max == 0);

    double fD_step = fs / (2.0 * N);
    int doppler_freq_size = static_cast<int>(fD_max / fD_step);
    int no_sub_tasks = N / r_max;

    RDMap mx(2 * doppler_freq_size + 1, IQBuffer(r_max, complex(0, 0)));

    // Reshape into sub-tasks
    mat ref_align = reshape(ref_ch, no_sub_tasks, r_max);
    mat surv_align = reshape(surv_ch, no_sub_tasks, r_max);

    for (auto& row : ref_align)
        row.resize(row.size() * 2);

    for (int i = 0; i < no_sub_tasks - 1; ++i) {
        int row_size = static_cast<int>(surv_align[i].size());
        surv_align[i].resize(row_size * 2);
        for (int j = row_size; j < row_size * 2; ++j)
            surv_align[i][j] = surv_align[i + 1][j - row_size];
    }
    surv_align.back().resize(surv_align.back().size() * 2, complex(0, 0));

    for (auto& row : ref_align)
        fft_backend_->forward(row);
    for (auto& row : surv_align)
        fft_backend_->forward(row);

    mat ref_conj = conjugate(ref_align);
    mat corr = multiply(surv_align, ref_conj);

    for (auto& row : corr)
        fft_backend_->inverse(row);

    // Transpose correlation matrix
    mat corr_t(2 * r_max, IQBuffer(2 * no_sub_tasks));
    assert(2 * r_max == static_cast<int>(corr[0].size()));
    assert(no_sub_tasks == static_cast<int>(corr.size()));
    for (int i = 0; i < no_sub_tasks; ++i)
        for (int j = 0; j < 2 * r_max; ++j)
            corr_t[j][i] = corr[i][j];

    for (auto& col : corr_t)
        fft_backend_->forward(col);

    for (int i = 0; i < doppler_freq_size; ++i) {
        for (int j = 0; j < r_max; ++j) {
            mx[i][j] = corr_t[j][2 * no_sub_tasks - doppler_freq_size + i];
            mx[doppler_freq_size + i][j] = corr_t[j][i];
        }
    }

    for (int j = 0; j < r_max; ++j)
        mx[doppler_freq_size][j] = corr_t[j][doppler_freq_size];

    return mx;
}

} // namespace hampr
