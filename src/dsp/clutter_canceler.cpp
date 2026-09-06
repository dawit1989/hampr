#include <hampr/dsp/clutter_canceler.hpp>
#include <hampr/core/exception.hpp>
#include <hampr/utils/fftw_wrapper.hpp>
#include <hampr/utils/linalg_utils.hpp>
#include <cmath>
#include <cassert>

namespace hampr {

IQBuffer ClutterCanceler::filter(const IQBuffer& ref_ch, const IQBuffer& surv_ch) const {
    return wiener_smi_mre(ref_ch, surv_ch);
}

IQBuffer ClutterCanceler::wiener_smi_mre(const IQBuffer& ref_ch, const IQBuffer& surv_ch) const {
    if (ref_ch.size() != surv_ch.size())
        throw HamprException("Mismatching channel sizes");

    int N = static_cast<int>(ref_ch.size());
    int K = K_;

    // Pruned autocorrelation of reference
    auto pruned_corr = [&](const IQBuffer& y, const IQBuffer& x, int clen) -> IQBuffer {
        int n = static_cast<int>(x.size());
        if (n != static_cast<int>(y.size()) || clen > (n + 1) / 2) {
                throw HamprException("Invalid sized arrays in pruned_corr");
        }

        int cols = clen - 1;
        int rows = n / cols + 1;
        int zeropad = cols * rows - n;

        IQBuffer x_pad = x;
        IQBuffer y_pad = y;
        for (int i = 0; i < zeropad; ++i) {
            x_pad.push_back(0);
            y_pad.push_back(0);
        }

        std::vector<IQBuffer> xp(rows), yp(rows);
        int cnt = 0;
        for (int i = 0; i < n; i += cols) {
            xp[cnt] = IQBuffer(x_pad.begin() + i, x_pad.begin() + i + cols);
            yp[cnt++] = IQBuffer(y_pad.begin() + i, y_pad.begin() + i + cols);
        }

        IQBuffer zeros(cols);
        for (int i = 0; i < rows - 1; ++i)
            yp[i].insert(yp[i].end(), yp[i + 1].begin(), yp[i + 1].end());
        yp.back().insert(yp.back().end(), zeros.begin(), zeros.end());
        for (int i = 0; i < rows; ++i)
            xp[i].insert(xp[i].end(), zeros.begin(), zeros.end());

        for (int i = 0; i < rows; ++i) {
            fft(xp[i]);
            fft(yp[i]);
        }

        auto bpw_conj = conjugate(yp);
        mat mult_result(rows, array(xp[0].size()));
        for (int i = 0; i < rows; ++i)
            for (size_t j = 0; j < xp[i].size(); ++j)
                mult_result[i][j] = xp[i][j] * bpw_conj[i][j];

        for (int i = 0; i < rows; ++i)
            ifft(mult_result[i]);

        mat fft_shifted = conjugate(fftshift(mult_result));

        IQBuffer output(clen);
        for (size_t i = 0; i < fft_shifted.size(); ++i) {
            int idx = 0;
            for (int j = clen - 1; j >= 0; --j)
                output[idx++] += fft_shifted[i][j];
        }
        return output;
    };

    IQBuffer pruned = pruned_corr(ref_ch, ref_ch, K);

    mat R(K, array(K));
    for (int i = 0; i < K; ++i)
        R[i][0] = pruned[i];

    IQBuffer r = pruned_corr(surv_ch, ref_ch, K);

    for (int k = 1; k < K; ++k) {
        IQBuffer shifted = shift(pruned, k);
        for (int i = 0; i < K; ++i)
            R[i][k] = shifted[i];
    }

    R = add(R, transposeConj(R));
    for (int i = 0; i < K; ++i)
        R[i][i] *= 0.5;

    IQBuffer w = dot(inverse(R), r);

    IQBuffer output = surv_ch;
    IQBuffer convolved = convolve(ref_ch, w);
    for (int i = 0; i < N; ++i)
        output[i] -= convolved[i];

    return output;
}

} // namespace hampr



