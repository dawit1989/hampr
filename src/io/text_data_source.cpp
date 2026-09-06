#include <hampr/io/text_data_source.hpp>
#include <fstream>
#include <hampr/core/exception.hpp>

namespace hampr {

IQMatrix TextDataSource::load(const std::string& filename, double& fs) {
    std::ifstream infile(filename);
    if (!infile.is_open())
        throw HamprException("Cannot open dataset file: " + filename);;

    infile >> fs;

    int rows, cols;
    infile >> rows >> cols;

    IQMatrix data(rows, IQBuffer(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) {
            double real, imag;
            infile >> real >> imag;
            data[i][j] = complex(real, imag);
        }

    return data;
}

TargetTrack TextDataSource::load_track(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open())
        throw HamprException("Cannot open target reference track: " + filename);;

    int rows, cols;
    infile >> rows >> cols;

    TargetTrack track;
    track.reserve(rows);
    for (int i = 0; i < rows; ++i) {
        TargetTrackPoint pt;
        double time_idx;
        infile >> time_idx >> pt.timestamp >> pt.latitude >> pt.longitude
               >> pt.altitude >> pt.speed >> pt.direction >> pt.range
               >> pt.doppler >> pt.azimuth;
        pt.time_index = static_cast<int>(time_idx);
        track.push_back(pt);
    }

    return track;
}

} // namespace hampr
