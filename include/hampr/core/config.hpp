#ifndef HAMPR_CORE_CONFIG_HPP
#define HAMPR_CORE_CONFIG_HPP

#include <string>
#include <vector>

namespace hampr {

struct Config {
    int ref_channel_index = 0;
    int filter_taps = 128;
    int num_antennas = 3;
    double antenna_spacing = 0.528;
    std::string window_type = "Hann";
    int search_window_size = 8;
    std::vector<int> metric_win = {6, 6, 3, 3};
    std::string td_filter_method = "wiener_smi_mre";
};

} // namespace hampr

#endif // HAMPR_CORE_CONFIG_HPP
