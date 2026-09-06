#ifndef HAMPR_DSP_DOA_ESTIMATOR_HPP
#define HAMPR_DSP_DOA_ESTIMATOR_HPP

#include <hampr/core/types.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <hampr/utils/math_utils.hpp>
#include <hampr/utils/linalg_utils.hpp>
#include <vector>

namespace hampr {

class DOAEstimator {
public:
    DOAEstimator() = default;

    std::vector<double> estimate(const std::vector<RDMap>& rd_maps,
                                 const std::vector<std::pair<int,int>>& hit_list,
                                 const std::vector<double>& array_alignment) const;

private:
    mat gen_ula_scanning_vectors(const std::vector<double>& alignment,
                                 const std::vector<double>& thetas) const;

    FlatBuffer gen_ula_scanning_vectors_fb(const std::vector<double>& alignment,
                                           const std::vector<double>& thetas) const;

    array doa_music(const mat& R,
                    const mat& scanning_vectors,
                    int signal_dimension) const;

    // Flat buffer variant for GPU acceleration
    array doa_music_fb(const FlatBuffer& R,
                       const FlatBuffer& scanning_vectors,
                       int signal_dimension) const;
};

} // namespace hampr

#endif // HAMPR_DSP_DOA_ESTIMATOR_HPP
