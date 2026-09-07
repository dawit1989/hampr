#include <hampr/dsp/array_geometry.hpp>
#include <hampr/core/exception.hpp>
#include <cmath>

namespace hampr {

ArrayGeometry::ArrayGeometry(const ArrayGeometrySpec& spec) {
    type_ = spec.type;
    element_spacing_ = spec.element_spacing;
    radius_ = spec.radius;

    switch (spec.type) {
        case ArrayGeometrySpec::Type::ULA: {
            int n = spec.ula_elements > 0 ? spec.ula_elements : 4;
            // Center the array: positions at -d*(n-1)/2, ..., 0, ..., d*(n-1)/2
            double d = spec.element_spacing;
            for (int i = 0; i < n; ++i) {
                double pos = d * (static_cast<double>(i) - (n - 1) / 2.0);
                elements_.push_back({pos, 0.0, 0.0});
            }
            break;
        }
        case ArrayGeometrySpec::Type::UCA: {
            int n = spec.ula_elements > 0 ? spec.ula_elements : 4;
            double r = spec.radius;
            for (int i = 0; i < n; ++i) {
                double phi = 2.0 * M_PI * static_cast<double>(i) / n;
                elements_.push_back({r * std::cos(phi), r * std::sin(phi), 0.0});
            }
            break;
        }
        case ArrayGeometrySpec::Type::URA: {
            int rows = spec.ura_rows > 0 ? spec.ura_rows : 2;
            int cols = spec.ura_cols > 0 ? spec.ura_cols : 2;
            double d = spec.element_spacing;
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c) {
                    double x = d * (static_cast<double>(c) - (cols - 1) / 2.0);
                    double y = d * (static_cast<double>(r) - (rows - 1) / 2.0);
                    elements_.push_back({x, y, 0.0});
                }
            break;
        }
        case ArrayGeometrySpec::Type::ARBITRARY: {
            for (const auto& el : spec.elements)
                elements_.push_back(el);
            break;
        }
    }

    if (elements_.empty())
        throw HamprException("ArrayGeometry: no elements configured");
}

mat ArrayGeometry::gen_scanning_vectors(const std::vector<double>& thetas,
                                        const std::vector<double>& phis) const {
    int M = static_cast<int>(elements_.size());
    int N = static_cast<int>(thetas.size());
    mat scanning_vectors(M, array(N));

    for (int i = 0; i < N; ++i) {
        double theta = thetas[i] * DEG_TO_RAD;
        double phi = (phis.size() == static_cast<size_t>(N))
                     ? phis[i] * DEG_TO_RAD : 0.0;
        for (int j = 0; j < M; ++j) {
            double phase = 0.0;
            switch (type_) {
                case ArrayGeometrySpec::Type::ULA:
                    phase = 2.0 * M_PI * elements_[j].x * std::cos(theta);
                    break;
                case ArrayGeometrySpec::Type::UCA:
                    phase = 2.0 * M_PI * elements_[j].x * std::cos(theta)
                          + 2.0 * M_PI * elements_[j].y * std::sin(theta);
                    break;
                case ArrayGeometrySpec::Type::URA:
                    phase = 2.0 * M_PI * (elements_[j].x * std::sin(theta) * std::cos(phi)
                                        + elements_[j].y * std::sin(theta) * std::sin(phi));
                    break;
                case ArrayGeometrySpec::Type::ARBITRARY:
                    phase = 2.0 * M_PI * (elements_[j].x * std::sin(theta) * std::cos(phi)
                                        + elements_[j].y * std::sin(theta) * std::sin(phi)
                                        + elements_[j].z * std::cos(theta));
                    break;
            }
            scanning_vectors[j][i] = complex(std::cos(phase), std::sin(phase));
        }
    }
    return scanning_vectors;
}

FlatBuffer ArrayGeometry::gen_scanning_vectors_fb(const std::vector<double>& thetas,
                                                  const std::vector<double>& phis) const {
    int M = static_cast<int>(elements_.size());
    int N = static_cast<int>(thetas.size());
    FlatBuffer scanning_vectors(M, N);

    for (int i = 0; i < N; ++i) {
        double theta = thetas[i] * DEG_TO_RAD;
        double phi = (phis.size() == static_cast<size_t>(N))
                     ? phis[i] * DEG_TO_RAD : 0.0;
        for (int j = 0; j < M; ++j) {
            double phase = 0.0;
            switch (type_) {
                case ArrayGeometrySpec::Type::ULA:
                    phase = 2.0 * M_PI * elements_[j].x * std::cos(theta);
                    break;
                case ArrayGeometrySpec::Type::UCA:
                    phase = 2.0 * M_PI * elements_[j].x * std::cos(theta)
                          + 2.0 * M_PI * elements_[j].y * std::sin(theta);
                    break;
                case ArrayGeometrySpec::Type::URA:
                    phase = 2.0 * M_PI * (elements_[j].x * std::sin(theta) * std::cos(phi)
                                        + elements_[j].y * std::sin(theta) * std::sin(phi));
                    break;
                case ArrayGeometrySpec::Type::ARBITRARY:
                    phase = 2.0 * M_PI * (elements_[j].x * std::sin(theta) * std::cos(phi)
                                        + elements_[j].y * std::sin(theta) * std::sin(phi)
                                        + elements_[j].z * std::cos(theta));
                    break;
            }
            scanning_vectors(j, i) = complex(std::cos(phase), std::sin(phase));
        }
    }
    return scanning_vectors;
}

array ArrayGeometry::manifold(double theta, double phi) const {
    double theta_rad = theta * DEG_TO_RAD;
    double phi_rad = phi * DEG_TO_RAD;
    int M = static_cast<int>(elements_.size());
    array result(M);
    for (int j = 0; j < M; ++j) {
        double phase = 0.0;
        switch (type_) {
            case ArrayGeometrySpec::Type::ULA:
                phase = 2.0 * M_PI * elements_[j].x * std::cos(theta_rad);
                break;
            case ArrayGeometrySpec::Type::UCA:
                phase = 2.0 * M_PI * elements_[j].x * std::cos(theta_rad)
                      + 2.0 * M_PI * elements_[j].y * std::sin(theta_rad);
                break;
            case ArrayGeometrySpec::Type::URA:
                phase = 2.0 * M_PI * (elements_[j].x * std::sin(theta_rad) * std::cos(phi_rad)
                                    + elements_[j].y * std::sin(theta_rad) * std::sin(phi_rad));
                break;
            case ArrayGeometrySpec::Type::ARBITRARY:
                phase = 2.0 * M_PI * (elements_[j].x * std::sin(theta_rad) * std::cos(phi_rad)
                                    + elements_[j].y * std::sin(theta_rad) * std::sin(phi_rad)
                                    + elements_[j].z * std::cos(theta_rad));
                break;
        }
        result[j] = complex(std::cos(phase), std::sin(phase));
    }
    return result;
}

} // namespace hampr

