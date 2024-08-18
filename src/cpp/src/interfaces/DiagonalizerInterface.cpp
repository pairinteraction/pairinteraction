#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>

namespace pairinteraction {
template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerInterface<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                    real_t min_eigenvalue, real_t max_eigenvalue,
                                    int precision) const {
    int dim = matrix.rows();

    auto eigensys = eigh(matrix, precision);

    auto *it_begin = std::lower_bound(eigensys.eigenvalues.data(),
                                      eigensys.eigenvalues.data() + dim, min_eigenvalue);
    auto *it_end = std::upper_bound(eigensys.eigenvalues.data(), eigensys.eigenvalues.data() + dim,
                                    max_eigenvalue);
    eigensys.eigenvectors = eigensys.eigenvectors
                                .block(0, std::distance(eigensys.eigenvalues.data(), it_begin), dim,
                                       std::distance(it_begin, it_end))
                                .eval();
    eigensys.eigenvalues = eigensys.eigenvalues
                               .segment(std::distance(eigensys.eigenvalues.data(), it_begin),
                                        std::distance(it_begin, it_end))
                               .eval();

    return eigensys;
}

// Explicit instantiations
template class DiagonalizerInterface<float>;
template class DiagonalizerInterface<double>;
template class DiagonalizerInterface<std::complex<float>>;
template class DiagonalizerInterface<std::complex<double>>;
} // namespace pairinteraction
