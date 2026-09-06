#include <hampr/dsp/channel_isolator.hpp>
#include <hampr/core/exception.hpp>
#include <iostream>

namespace hampr {

std::pair<IQBuffer, IQMatrix>
ChannelIsolator::isolate(const IQMatrix& iq_samples) const {
    int N = static_cast<int>(iq_samples[0].size());
    int Mp1 = static_cast<int>(iq_samples.size());

    if (Mp1 > N) {
        std::cerr << "WARNING: Number of samples exceeds antenna channels; input may be flipped\n";
    }

    if (ref_index_ < 0 || ref_index_ >= Mp1)
        throw HamprException("Invalid ref_index");

    IQBuffer ref_channel = iq_samples[ref_index_];
    IQMatrix surv_channels(Mp1 - 1);
    int idx = 0;
    for (int i = 0; i < Mp1; ++i)
        if (i != ref_index_)
            surv_channels[idx++] = iq_samples[i];

    return {ref_channel, surv_channels};
}

} // namespace hampr
