#ifndef HAMPR_DSP_CHANNEL_ISOLATOR_HPP
#define HAMPR_DSP_CHANNEL_ISOLATOR_HPP

#include <hampr/core/types.hpp>
#include <hampr/utils/math_utils.hpp>
#include <utility>

namespace hampr {

class ChannelIsolator {
public:
    ChannelIsolator(int ref_channel_index = 0)
        : ref_index_(ref_channel_index) {}

    std::pair<IQBuffer, IQMatrix>
    isolate(const IQMatrix& iq_samples) const;

private:
    int ref_index_;
};

} // namespace hampr

#endif // HAMPR_DSP_CHANNEL_ISOLATOR_HPP
