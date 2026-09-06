#include <hampr/dsp/doa_estimator.hpp>
#include <hampr/accel/linalg_backend.hpp>
#include <hampr/core/exception.hpp>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <algorithm>

namespace hampr {

mat DOAEstimator::gen_ula_scanning_vectors(
        const std::vector<double>& alignment,
        const std::vector<double>& thetas) const {
    int thetas_size = static_cast<int>(thetas.size());
    int M = static_cast<int>(alignment.size());
    mat scanning_vectors(M, array(thetas_size));
    for (int i = 0; i < thetas_size; ++i)
        for (int j = 0; j < M; ++j) {
            double phase = 2.0 * M_PI * alignment[j] * std::cos(0.0174533 * thetas[i]);
            scanning_vectors[j][i] = complex(std::cos(phase), std::sin(phase));
        }
    return scanning_vectors;
}

FlatBuffer DOAEstimator::gen_ula_scanning_vectors_fb(
        const std::vector<double>& alignment,
        const std::vector<double>& thetas) const {
    int thetas_size = static_cast<int>(thetas.size());
    int M = static_cast<int>(alignment.size());
    FlatBuffer scanning_vectors(M, thetas_size);
    for (int i = 0; i < thetas_size; ++i)
        for (int j = 0; j < M; ++j) {
            double phase = 2.0 * M_PI * alignment[j] * std::cos(0.0174533 * thetas[i]);
            scanning_vectors(j, i) = complex(std::cos(phase), std::sin(phase));
        }
    return scanning_vectors;
}

array DOAEstimator::doa_music(const mat& R,
                              const mat& scanning_vectors,
                              int signal_dimension) const {
    if (R.empty() || R.size() != R[0].size())
        throw HamprException("Correlation matrix is not square");
    if (R.size() != scanning_vectors.size())
        throw HamprException("Dimension mismatch between R and scanning_vectors");

    int R_size = static_cast<int>(R.size());
    Eigen::MatrixXcd mat_R(R_size, R_size);
    for (int i = 0; i < R_size; ++i)
        for (int j = 0; j < R_size; ++j)
            mat_R.row(i)[j] = R[i][j];

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(mat_R);
    Eigen::VectorXcd vals = ces.eigenvalues();
    Eigen::MatrixXcd vecs = ces.eigenvectors();

    int M = R_size;
    int noise_dimension = M - signal_dimension;

    // Sort by eigenvalue magnitude (ascending)
    std::vector<std::pair<double, array>> eig_pairs(M);
    for (int i = 0; i < M; ++i) {
        eig_pairs[i].first = vals[i].real() + vals[i].imag();
        eig_pairs[i].second = array(M);
        for (int j = 0; j < M; ++j)
            eig_pairs[i].second[j] = vecs.row(j)[i];
    }
    std::sort(eig_pairs.begin(), eig_pairs.end(),
              [](const auto& a, const auto& b) {
                  return a.first < b.first;
              });

    // Noise subspace (smallest eigenvalues)
    mat E(M, array(noise_dimension));
    for (int i = 0; i < noise_dimension; ++i)
        for (int j = 0; j < M; ++j)
            E[j][i] = eig_pairs[i].second[j];

    int scanning_cols = static_cast<int>(scanning_vectors[0].size());
    array ADORT(scanning_cols, 0);

    // Hoist the constant noise subspace projector E * E^H(out of the scan loop
    mat E_Et = matMult(E, transposeConj(E));

    for (int i = 0; i < scanning_cols; ++i) {
        complex denom(0.0, 0.0);
        for (int j = 0; j < M; ++j) {
            complex temp(0.0, 0.0);
            for (int k = 0; k < M; ++k)
                temp += E_Et[j][k] * scanning_vectors[k][i];
            denom += std::conj(scanning_vectors[j][i]) * temp;
        }
        ADORT[i] = 1.0 / std::abs(denom);
    }

    return ADORT;
}

array DOAEstimator::doa_music_fb(const FlatBuffer& R,
                                 const FlatBuffer& scanning_vectors,
                                 int signal_dimension) const {
    auto backend = create_linalg_backend();

    // Eigen decomposition via backend flat buffer method
    auto [eigenvalues, eigenvectors_fb] = backend->eigen_decompose_fb(R);
    mat eigenvectors = eigenvectors_fb.to_mat();

    int M = static_cast<int>(R.rows());
    int noise_dimension = M - signal_dimension;

    // Sort by eigenvalue magnitude (ascending)
    std::vector<std::pair<double, array>> eig_pairs(M);
    for (int i = 0; i < M; ++i) {
        eig_pairs[i].first = eigenvalues[i].real() + eigenvalues[i].imag();
        eig_pairs[i].second = array(M);
        for (int j = 0; j < M; ++j)
            eig_pairs[i].second[j] = eigenvectors[j][i];
    }
    std::sort(eig_pairs.begin(), eig_pairs.end(),
              [](const auto& a, const auto& b) {
                  return a.first < b.first;
              });

    // Noise subspace (smallest eigenvalues)
    mat E(M, array(noise_dimension));
    for (int i = 0; i < noise_dimension; ++i)
        for (int j = 0; j < M; ++j)
            E[j][i] = eig_pairs[i].second[j];

    int scanning_cols = scanning_vectors.cols();
    array ADORT(scanning_cols, 0);

    // Hoist E * E^H out of the scan loop
    mat E_Et = matMult(E, transposeConj(E));

    // Convert scanning_vectors to mat for the scan loop
    mat scan_mat = scanning_vectors.to_mat();

    for (int i = 0; i < scanning_cols; ++i) {
        complex denom(0.0, 0.0);
        for (int j = 0; j < M; ++j) {
            complex temp(0.0, 0.0);
            for (int k = 0; k < M; ++k)
                temp += E_Et[j][k] * scan_mat[k][i];
            denom += std::conj(scan_mat[j][i]) * temp;
        }
        ADORT[i] = 1.0 / std::abs(denom);
    }

    return ADORT;
}

std::vector<double> DOAEstimator::estimate(
        const std::vector<RDMap>& rd_maps,
        const std::vector<std::pair<int,int>>& hit_list,
        const std::vector<double>& array_alignment) const {

    std::vector<double> thetas(180);
    for (int i = 0; i < 180; ++i)
        thetas[i] = static_cast<double>(i);

    mat scanning_vectors = gen_ula_scanning_vectors(array_alignment, thetas);
    FlatBuffer scanning_vectors_fb = gen_ula_scanning_vectors_fb(array_alignment, thetas);

    std::vector<double> doa_list;
    for (const auto& hit : hit_list) {
        array az_vector(rd_maps.size());
        for (size_t j = 0; j < rd_maps.size(); ++j)
            az_vector[j] = rd_maps[j][hit.second][hit.first];

        array az_conj = conjugate(az_vector);
        mat R(az_conj.size(), array(az_conj.size()));
        for (size_t j = 0; j < az_conj.size(); ++j)
            for (size_t k = 0; k < az_conj.size(); ++k)
                R[j][k] = az_vector[j] * az_conj[k];

        // Use flat buffer variant for GPU-ready computation path
        FlatBuffer R_fb;
        R_fb.from_mat(R);

        array doa_res = doa_music_fb(R_fb, scanning_vectors_fb, 1);

        double max_val = 0;
        int argmax = 0;
        for (size_t j = 0; j < doa_res.size(); ++j) {
            double val = doa_res[j].real() + doa_res[j].imag();
            if (val > max_val) {
                max_val = val;
                argmax = static_cast<int>(j);
            }
        }
        doa_list.push_back(static_cast<double>(argmax));
    }
    return doa_list;
}

} // namespace hampr
