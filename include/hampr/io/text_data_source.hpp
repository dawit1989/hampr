#ifndef HAMPR_IO_TEXT_DATA_SOURCE_HPP
#define HAMPR_IO_TEXT_DATA_SOURCE_HPP

#include <hampr/io/data_source.hpp>

namespace hampr {

class TextDataSource : public DataSource {
public:
    IQMatrix load(const std::string& filename, double& fs) override;
    TargetTrack load_track(const std::string& filename) override;
};

} // namespace hampr

#endif // HAMPR_IO_TEXT_DATA_SOURCE_HPP
