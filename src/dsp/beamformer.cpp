#include <hampr/dsp/beamformer.hpp>
#include <hampr/core/exception.hpp>
#include <hampr/utils/linalg_utils.hpp>
#include <cmath>

namespace hampr {

Beamformer::Result Beamformer::beamform(
        const IQMatrix& ant_channels,
        const std::vector<double>& alignment,
        double direction) const {

    int N = static_cast<int>(ant_channels.size());
    int M = static_cast<int>(ant_channels[0].size());
    int align_size = static_cast<int>(alignment.size());

    if (M > N)
        throw HamprException("M > N: antenna channels exceed samples");
    if (M != align_size)
        throw HamprException("Mismatch in channel numbers");

    // Estimate spatial correlation matrix: R = X^T * X^H / N
    mat R(M, array(M));
    mat t = transpose(ant_channels);
    R = matMult(t, transposeConj(t));
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            R[i][j] /= complex(static_cast<double>(N), 0);

    // Array response vector: aS[j] = exp(1j * 2*pi * d * cos(theta))
    mat aS(align_size, array(1));
    double theta = 0.0174533 * direction;
    for (int j = 0; j < align_size; ++j) {
        double phase = 2.0 * M_PI * alignment[j] * std::cos(theta);
        aS[j][0] = complex(std::cos(phase), std::sin(phase));
    }

    // Wiener beamformer: w = R^{-1} * aS
    mat w = matMult(inverse(R), aS);

    // Normalize
    double norm_sum = 0;
    for (size_t j = 0; j < w.size(); ++j)
        for (size_t k = 0; k < w[j].size(); ++k)
            norm_sum += std::abs(w[j][k]);
    double norm_factor = std::sqrt(static_cast<double>(M)) / norm_sum;
    for (size_t j = 0; j < w.size(); ++j)
        for (size_t k = 0; k < w[j].size(); ++k)
            w[j][k] *= norm_factor;

    // Beamform: surv = w^H * X^T
    mat surv = matMult(transposeConj(w), transpose(ant_channels));

    Result result;
    result.surveillance_channel = surv[0];
    result.weights.resize(w.size());
    for (size_t j = 0; j < w.size(); ++j)
        result.weights[j] = w[j][0];
    return result;
}

} // namespace hampr

