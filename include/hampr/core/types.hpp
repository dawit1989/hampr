#ifndef HAMPR_CORE_TYPES_HPP
#define HAMPR_CORE_TYPES_HPP

#include <complex>
#include <vector>

namespace hampr {

using complex = std::complex<double>;
using array = std::vector<complex>;
using mat = std::vector<array>;
using IQBuffer = std::vector<complex>;
using IQMatrix = std::vector<IQBuffer>;
using RDMap = std::vector<std::vector<complex>>;

struct TargetTrackPoint {
    int time_index;
    double timestamp;
    double latitude;
    double longitude;
    double altitude;
    double speed;
    double direction;
    double range;
    double doppler;
    double azimuth;
};

using TargetTrack = std::vector<TargetTrackPoint>;

struct ProcessingResult {
    int target_range_cell_ref;
    int target_range_cell_found;
    int target_doppler_cell_ref;
    int target_doppler_cell_found;
    double azimuth_ref;
    double azimuth_estimated;
    double sinr;
    double isolation_time;
    double beamform_time;
    double filter_time;
    double detector_time;
    double target_find_time;
    double doa_time;
    double total_time;
};

} // namespace hampr

#endif // HAMPR_CORE_TYPES_HPP
