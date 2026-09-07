#ifndef HAMPR_DSP_ARRAY_GEOMETRY_HPP
#define HAMPR_DSP_ARRAY_GEOMETRY_HPP

#include <hampr/core/multi_site_types.hpp>
#include <hampr/core/flat_buffer.hpp>
#include <vector>

namespace hampr {

class ArrayGeometry {
public:
    ArrayGeometry() = default;
    explicit ArrayGeometry(const ArrayGeometrySpec& spec);

    int num_elements() const { return static_cast<int>(elements_.size()); }
    const std::vector<AntennaElement>& elements() const { return elements_; }
    ArrayGeometrySpec::Type type() const { return type_; }

    mat gen_scanning_vectors(const std::vector<double>& thetas,
                             const std::vector<double>& phis = {}) const;

    FlatBuffer gen_scanning_vectors_fb(const std::vector<double>& thetas,
                                       const std::vector<double>& phis = {}) const;

    array manifold(double theta, double phi = 0.0) const;

private:
    std::vector<AntennaElement> elements_;
    ArrayGeometrySpec::Type type_ = ArrayGeometrySpec::Type::ULA;
    double element_spacing_ = 0.5;
    double radius_ = 0.0;

    static constexpr double DEG_TO_RAD = 0.017453292519943295;
};

} // namespace hampr

#endif // HAMPR_DSP_ARRAY_GEOMETRY_HPP
