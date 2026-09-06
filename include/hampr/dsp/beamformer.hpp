#ifndef HAMPR_DSP_BEAMFORMER_HPP
#define HAMPR_DSP_BEAMFORMER_HPP

#include <hampr/core/types.hpp>
#include <hampr/utils/math_utils.hpp>
#include <vector>

namespace hampr {

class Beamformer {
public:
    struct Result {
        IQBuffer surveillance_channel;
        std::vector<complex> weights;
    };

    Beamformer(double antenna_spacing = 0.528)
        : spacing_(antenna_spacing) {}

    Result beamform(const IQMatrix& ant_channels,
                    const std::vector<double>& alignment,
                    double direction) const;

private:
    double spacing_;
};

} // namespace hampr

#endif // HAMPR_DSP_BEAMFORMER_HPP
