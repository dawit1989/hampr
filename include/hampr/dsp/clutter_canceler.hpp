#ifndef HAMPR_DSP_CLUTTER_CANCELER_HPP
#define HAMPR_DSP_CLUTTER_CANCELER_HPP

#include <hampr/core/types.hpp>
#include <hampr/utils/math_utils.hpp>
#include <string>

namespace hampr {

class ClutterCanceler {
public:
    explicit ClutterCanceler(int filter_taps = 128)
        : K_(filter_taps) {}

    IQBuffer filter(const IQBuffer& ref_ch,
                    const IQBuffer& surv_ch) const;

private:
    int K_;
    IQBuffer wiener_smi_mre(const IQBuffer& ref_ch,
                            const IQBuffer& surv_ch) const;
};

} // namespace hampr

#endif // HAMPR_DSP_CLUTTER_CANCELER_HPP
