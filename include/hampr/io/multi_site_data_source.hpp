#ifndef HAMPR_IO_MULTI_SITE_DATA_SOURCE_HPP
#define HAMPR_IO_MULTI_SITE_DATA_SOURCE_HPP

#include <hampr/io/data_source.hpp>
#include <hampr/io/streaming_source.hpp>
#include <hampr/core/multi_site_types.hpp>
#include <string>
#include <map>
#include <memory>

namespace hampr {

class MultiSiteDataSource : public DataSource {
public:
    MultiSiteDataSource() = default;
    ~MultiSiteDataSource() override = default;

    bool load_site(const std::string& site_id, const std::string& filename);
    bool load_metadata(const std::string& metadata_file);

    MultiSiteData get_multi_site_data();

    bool open_site(const std::string& site_id, const std::string& filename);
    IQMatrix next_batch_site(const std::string& site_id, size_t batch_size);
    bool has_more_site(const std::string& site_id) const;
    double sampling_rate_site(const std::string& site_id) const;
    int num_channels_site(const std::string& site_id) const;
    void reset_all();
    void close_all();

    IQMatrix load(const std::string& filename, double& fs) override;
    TargetTrack load_track(const std::string& filename) override;

    const std::vector<ReceiverInfo>& receivers() const { return receiver_infos_; }
    double sampling_rate() const { return master_fs_; }
    bool empty() const { return site_iq_data_.empty(); }

private:
    std::map<std::string, IQMatrix> site_iq_data_;
    std::map<std::string, std::unique_ptr<StreamingDataSource>> site_sources_;
    std::vector<ReceiverInfo> receiver_infos_;
    double master_fs_ = 0.0;

    static std::string find_file(const std::string& filename);
};

} // namespace hampr

#endif // HAMPR_IO_MULTI_SITE_DATA_SOURCE_HPP
