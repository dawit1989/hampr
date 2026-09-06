#include <hampr/dsp/metric_extractor.hpp>
#include <cmath>

namespace hampr {

double MetricExtractor::extract_snr(RDMap& rd_matrix,
                                    std::vector<int> target_rd,
                                    std::vector<int>& win) const {
    target_rd[1] += (static_cast<int>(rd_matrix.size()) - 1) / 2;

    mat mask(2 * win[1] + 1, IQBuffer(2 * win[0] + 1, complex(1, 0)));
    double cell_counter = static_cast<double>(mask.size() * mask[0].size());

    for (int i = win[1] - win[3]; i < win[1] + win[3] + 1; ++i)
        for (int j = win[0] - win[2]; j < win[0] + win[2] + 1; ++j) {
            mask[i][j] = complex(0, 0);
            --cell_counter;
        }

    mat rd_block(2 * win[1] + 1);
    for (int i = 0; i < 2 * win[1] + 1; ++i)
        for (int j = 0; j < 2 * win[0] + 1; ++j)
            rd_block[i].push_back(rd_matrix[i + target_rd[1] - win[1]][j + target_rd[0] - win[0]]);

    rd_block = multiply(rd_block, mask);
    rd_block = multiply(rd_block, mask);

    for (auto& row : rd_block)
        for (auto& col : row)
            col = complex(std::abs(col) * std::abs(col), 0);

    double P_env = sum(rd_block).real() / cell_counter;
    double P_target = std::abs(rd_matrix[target_rd[1]][target_rd[0]]) *
                      std::abs(rd_matrix[target_rd[1]][target_rd[0]]);

    return 10.0 * std::log10(P_target / P_env);
}

} // namespace hampr
