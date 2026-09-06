#ifndef HAMPR_IO_DATA_SOURCE_HPP
#define HAMPR_IO_DATA_SOURCE_HPP

#include <hampr/core/types.hpp>
#include <string>

namespace hampr {

class DataSource {
public:
    virtual ~DataSource() = default;

    virtual IQMatrix load(const std::string& filename, double& fs) = 0;
    virtual TargetTrack load_track(const std::string& filename) = 0;
};

} // namespace hampr

#endif // HAMPR_IO_DATA_SOURCE_HPP
