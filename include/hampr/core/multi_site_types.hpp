#ifndef HAMPR_CORE_MULTI_SITE_TYPES_HPP
#define HAMPR_CORE_MULTI_SITE_TYPES_HPP

#include <hampr/core/types.hpp>
#include <string>
#include <vector>

namespace hampr {

struct AntennaElement {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct ArrayGeometrySpec {
    enum class Type { ULA, UCA, URA, ARBITRARY };
    Type type = Type::ULA;
    double element_spacing = 0.5;
    double radius = 0.0;
    int ula_elements = 4;
    int ura_rows = 0;
    int ura_cols = 0;
    std::vector<AntennaElement> elements;
};

struct ReceiverInfo {
    std::string id;
    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double clock_offset = 0.0;
    double gain = 0.0;
    double phase_offset = 0.0;
    ArrayGeometrySpec array_geometry;
    int ref_channel_index = 0;
    int num_channels = 4;
};

struct MultiSiteData {
    std::vector<ReceiverInfo> receivers;
    std::vector<IQMatrix> iq_data;
    double fs = 0.0;
    double master_timestamp = 0.0;
};

struct SynchronizedBatch {
    std::vector<ReceiverInfo> receivers;
    std::vector<IQMatrix> aligned_data;
    double fs = 0.0;
    int batch_samples = 0;
};

struct MultiSiteResult {
    std::vector<ProcessingResult> site_results;
    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double geolocation_error = 0.0;
    double fusion_time = 0.0;
    double ms_confidence = 0.0;
    double total_time = 0.0;
};

} // namespace hampr

#endif // HAMPR_CORE_MULTI_SITE_TYPES_HPP
