#include <hampr/io/multi_site_data_source.hpp>
#include <hampr/io/text_data_source.hpp>
#include <hampr/utils/json.hpp>
#include <hampr/core/exception.hpp>
#include <fstream>
#include <sstream>
#include <filesystem>

namespace hampr {

namespace fs = std::filesystem;

static std::string read_file(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
        throw HamprException("Cannot open file: " + filename);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::string MultiSiteDataSource::find_file(const std::string& filename) {
    if (fs::exists(filename))
        return filename;
    std::vector<std::string> candidates = {
        "dataset/" + filename,
        "../passiveradar_data/dataset/" + filename,
        "../../passiveradar_data/dataset/" + filename,
    };
    for (const auto& c : candidates) {
        if (fs::exists(c))
            return c;
    }
    return filename;
}

bool MultiSiteDataSource::load_site(const std::string& site_id, const std::string& filename) {
    try {
        std::string path = find_file(filename);
        TextDataSource loader;
        double fs = 0.0;
        IQMatrix data = loader.load(path, fs);
        site_iq_data_[site_id] = std::move(data);
        if (master_fs_ == 0.0)
            master_fs_ = fs;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool MultiSiteDataSource::load_metadata(const std::string& metadata_file) {
    try {
        std::string path = find_file(metadata_file);
        std::string content = read_file(path);
        auto root = json::parse(content);

        if (root.is_number())
            master_fs_ = root.as_number();
        else if (root.is_object()) {
            if (root.obj.count("fs"))
                master_fs_ = root.at("fs").as_number();
            if (root.obj.count("receivers")) {
                receiver_infos_.clear();
                for (size_t i = 0; i < root.at("receivers").arr.size(); ++i) {
                    const auto& r = root.at("receivers").at(i);
                    ReceiverInfo info;
                    if (r.obj.count("id"))        info.id = r.at("id").as_string();
                    if (r.obj.count("latitude"))   info.latitude = r.at("latitude").as_number();
                    if (r.obj.count("longitude"))  info.longitude = r.at("longitude").as_number();
                    if (r.obj.count("altitude"))   info.altitude = r.at("altitude").as_number();
                    if (r.obj.count("clock_offset")) info.clock_offset = r.at("clock_offset").as_number();
                    if (r.obj.count("gain"))       info.gain = r.at("gain").as_number();
                    if (r.obj.count("phase_offset")) info.phase_offset = r.at("phase_offset").as_number();
                    if (r.obj.count("ref_channel_index")) info.ref_channel_index = r.at("ref_channel_index").as_int();
                    if (r.obj.count("num_channels")) info.num_channels = r.at("num_channels").as_int();

                    if (r.obj.count("array_geometry")) {
                        const auto& ag = r.at("array_geometry");
                        if (ag.obj.count("type")) {
                            std::string t = ag.at("type").as_string();
                            if (t == "ULA") info.array_geometry.type = ArrayGeometrySpec::Type::ULA;
                            else if (t == "UCA") info.array_geometry.type = ArrayGeometrySpec::Type::UCA;
                            else if (t == "URA") info.array_geometry.type = ArrayGeometrySpec::Type::URA;
                            else info.array_geometry.type = ArrayGeometrySpec::Type::ARBITRARY;
                        }
                        if (ag.obj.count("element_spacing")) info.array_geometry.element_spacing = ag.at("element_spacing").as_number();
                        if (ag.obj.count("radius")) info.array_geometry.radius = ag.at("radius").as_number();
                        if (ag.obj.count("ula_elements")) info.array_geometry.ula_elements = ag.at("ula_elements").as_int();
                        if (ag.obj.count("ura_rows")) info.array_geometry.ura_rows = ag.at("ura_rows").as_int();
                        if (ag.obj.count("ura_cols")) info.array_geometry.ura_cols = ag.at("ura_cols").as_int();
                        if (ag.obj.count("elements") && ag.at("elements").is_array()) {
                            info.array_geometry.elements.clear();
                            for (size_t e = 0; e < ag.at("elements").arr.size(); ++e) {
                                const auto& el = ag.at("elements").at(e);
                                AntennaElement ae;
                                if (el.obj.count("x")) ae.x = el.at("x").as_number();
                                if (el.obj.count("y")) ae.y = el.at("y").as_number();
                                if (el.obj.count("z")) ae.z = el.at("z").as_number();
                                info.array_geometry.elements.push_back(ae);
                            }
                        }
                    }
                    receiver_infos_.push_back(info);
                }
            }
        }
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

MultiSiteData MultiSiteDataSource::get_multi_site_data() {
    MultiSiteData data;
    data.fs = master_fs_;

    int site_idx = 0;
    for (const auto& [id, iq] : site_iq_data_) {
        if (site_idx < static_cast<int>(receiver_infos_.size()))
            data.receivers.push_back(receiver_infos_[site_idx]);
        else {
            ReceiverInfo info;
            info.id = id;
            info.num_channels = static_cast<int>(iq.size());
            data.receivers.push_back(info);
        }
        data.iq_data.push_back(iq);
        ++site_idx;
    }
    return data;
}

bool MultiSiteDataSource::open_site(const std::string& site_id, const std::string& filename) {
    try {
        std::string path = find_file(filename);
        auto source = std::make_unique<StreamingDataSource>();
        if (!source->open(path)) {
            return false;
        }
        if (master_fs_ == 0.0)
            master_fs_ = source->sampling_rate();
        site_sources_[site_id] = std::move(source);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

IQMatrix MultiSiteDataSource::next_batch_site(const std::string& site_id, size_t batch_size) {
    auto it = site_sources_.find(site_id);
    if (it == site_sources_.end())
        throw HamprException("MultiSiteDataSource: unknown site: " + site_id);
    return it->second->next_batch(batch_size);
}

bool MultiSiteDataSource::has_more_site(const std::string& site_id) const {
    auto it = site_sources_.find(site_id);
    if (it == site_sources_.end())
        return false;
    return it->second->has_more();
}

double MultiSiteDataSource::sampling_rate_site(const std::string& site_id) const {
    auto it = site_sources_.find(site_id);
    if (it == site_sources_.end())
        return 0.0;
    return it->second->sampling_rate();
}

int MultiSiteDataSource::num_channels_site(const std::string& site_id) const {
    auto it = site_sources_.find(site_id);
    if (it == site_sources_.end())
        return 0;
    return it->second->num_channels();
}

void MultiSiteDataSource::reset_all() {
    for (auto& [id, source] : site_sources_)
        source->reset();
}

void MultiSiteDataSource::close_all() {
    for (auto& [id, source] : site_sources_)
        source->close();
    site_sources_.clear();
}

IQMatrix MultiSiteDataSource::load(const std::string& filename, double& fs) {
    TextDataSource loader;
    return loader.load(filename, fs);
}

TargetTrack MultiSiteDataSource::load_track(const std::string& filename) {
    TextDataSource loader;
    return loader.load_track(filename);
}

} // namespace hampr

