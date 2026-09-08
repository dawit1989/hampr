#ifndef HAMPR_DSP_CLUTTER_CANCELER_HPP
#define HAMPR_DSP_CLUTTER_CANCELER_HPP

#include <hampr/core/types.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <hampr/accel/fft_backend.hpp>
#include <hampr/accel/linalg_backend.hpp>
#include <hampr/utils/math_utils.hpp>
#include <string>
#include <memory>

namespace hampr {

class ClutterCanceler {
public:
    explicit ClutterCanceler(int filter_taps = 128)
        : K_(filter_taps),
          fft_backend_(create_fft_backend()),
          linalg_backend_(create_linalg_backend()) {}

    IQBuffer filter(const IQBuffer& ref_ch,
                    const IQBuffer& surv_ch) const;

private:
    int K_;
    std::unique_ptr<FFTBackend> fft_backend_;
    std::unique_ptr<LinAlgBackend> linalg_backend_;
    IQBuffer wiener_smi_mre(const IQBuffer& ref_ch,
                            const IQBuffer& surv_ch) const;
};

} // namespace hampr

#endif // HAMPR_DSP_CLUTTER_CANCELER_HPP
