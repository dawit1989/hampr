#ifndef HAMPR_DSP_DETECTOR_HPP
#define HAMPR_DSP_DETECTOR_HPP

#include <hampr/core/types.hpp>
#include <hampr/utils/math_utils.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <hampr/accel/fft_backend.hpp>
#include <hampr/accel/linalg_backend.hpp>
#include <hampr/utils/math_utils.hpp>
#include <string>
#include <vector>
#include <memory>

namespace hampr {

class Detector {
public:
    explicit Detector(const std::string& window_type = "Hann")
        : window_type_(window_type),
          fft_backend_(create_fft_backend()),
          linalg_backend_(create_linalg_backend()) {}

    IQBuffer apply_window(const IQBuffer& surv_ch) const;

    RDMap detect(IQBuffer& ref_ch,
                 IQBuffer& surv_ch,
                 double fs,
                 int fD_max,
                 int r_max) const;

    std::vector<double> hann(int M) const;

private:
    std::string window_type_;
    std::unique_ptr<FFTBackend> fft_backend_;
    std::unique_ptr<LinAlgBackend> linalg_backend_;

    RDMap cc_detector_ons(IQBuffer& ref_ch,
                          IQBuffer& surv_ch,
                          double fs,
                          int fD_max,
                          int r_max) const;
};

} // namespace hampr

#endif // HAMPR_DSP_DETECTOR_HPP
