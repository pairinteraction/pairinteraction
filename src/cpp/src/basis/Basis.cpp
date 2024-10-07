#include "pairinteraction/basis/Basis.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisClassicalLight.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetClassicalLight.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/wigner.hpp"

#include <numeric>
#include <set>

namespace pairinteraction {

// Forward declarations
template <typename Scalar>
class BasisAtom;
template <typename Scalar>
class BasisClassicalLight;

template <typename Derived>
void Basis<Derived>::perform_sorter_checks(const std::vector<TransformationType> &labels) const {
    // Check if the labels are valid sorting labels
    for (const auto &label : labels) {
        if (!utils::is_sorting(label)) {
            throw std::invalid_argument("One of the labels is not a valid sorting label.");
        }
    }
}

template <typename Derived>
void Basis<Derived>::perform_blocks_checks(
    const std::set<TransformationType> &unique_labels) const {
    // Check if the states are sorted by the requested labels
    std::set<TransformationType> unique_labels_present;
    for (const auto &label : get_transformation().transformation_type) {
        if (!utils::is_sorting(label) || unique_labels_present.size() >= unique_labels.size()) {
            break;
        }
        unique_labels_present.insert(label);
    }
    if (unique_labels != unique_labels_present) {
        throw std::invalid_argument("The states are not sorted by the requested labels.");
    }

    // Throw a meaningful error if getting the blocks by energy is requested as this might be a
    // common mistake
    if (unique_labels.count(TransformationType::SORT_BY_ENERGY) > 0) {
        throw std::invalid_argument("States do not store the energy and thus no energy blocks can "
                                    "be obtained. Use an energy operator instead.");
    }
}

template <typename Derived>
Basis<Derived>::Basis(ketvec_t &&kets)
    : kets(std::move(kets)), coefficients{{static_cast<Eigen::Index>(this->kets.size()),
                                           static_cast<Eigen::Index>(this->kets.size())},
                                          {TransformationType::SORT_BY_KET}} {
    if (this->kets.empty()) {
        throw std::invalid_argument("The basis must contain at least one element.");
    }
    quantum_number_f_of_states.reserve(this->kets.size());
    quantum_number_m_of_states.reserve(this->kets.size());
    parity_of_states.reserve(this->kets.size());
    ket_id_to_index.reserve(this->kets.size());
    size_t index = 0;
    for (const auto &ket : this->kets) {
        quantum_number_f_of_states.push_back(ket->get_quantum_number_f());
        quantum_number_m_of_states.push_back(ket->get_quantum_number_m());
        parity_of_states.push_back(ket->get_parity());
        ket_id_to_index[ket->get_id()] = index++;
    }
    ket_of_states.resize(this->kets.size());
    std::iota(ket_of_states.begin(), ket_of_states.end(), 0);
    coefficients.matrix.setIdentity();
}

template <typename Derived>
const Derived &Basis<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
const typename Basis<Derived>::ketvec_t &Basis<Derived>::get_kets() const {
    return kets;
}

template <typename Derived>
const Eigen::SparseMatrix<typename Basis<Derived>::scalar_t, Eigen::RowMajor> &
Basis<Derived>::get_coefficients() const {
    return coefficients.matrix;
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_quantum_number_f(size_t index_state) const {
    if (index_state >= static_cast<size_t>(coefficients.matrix.cols())) {
        throw std::out_of_range("The index is out of bounds.");
    }
    if (quantum_number_f_of_states[index_state] == std::numeric_limits<real_t>::max()) {
        throw std::invalid_argument("The state does not have a well-defined quantum number f.");
    }
    return quantum_number_f_of_states[index_state];
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_quantum_number_m(size_t index_state) const {
    if (index_state >= static_cast<size_t>(coefficients.matrix.cols())) {
        throw std::out_of_range("The index is out of bounds.");
    }
    if (quantum_number_m_of_states[index_state] == std::numeric_limits<real_t>::max()) {
        throw std::invalid_argument("The state does not have a well-defined quantum number m.");
    }
    return quantum_number_m_of_states[index_state];
}

template <typename Derived>
Parity Basis<Derived>::get_parity(size_t index_state) const {
    if (index_state >= static_cast<size_t>(coefficients.matrix.cols())) {
        throw std::out_of_range("The index is out of bounds.");
    }
    if (parity_of_states[index_state] == Parity::UNKNOWN) {
        throw std::invalid_argument("The state does not have a well-defined parity.");
    }
    return parity_of_states[index_state];
}

template <typename Derived>
std::shared_ptr<const typename Basis<Derived>::ket_t>
Basis<Derived>::get_ket_with_largest_overlap(size_t index_state) const {
    if (index_state >= static_cast<size_t>(coefficients.matrix.cols())) {
        throw std::out_of_range("The index is out of bounds.");
    }
    if (ket_of_states[index_state] == std::numeric_limits<int>::max()) {
        throw std::invalid_argument("The state does not belong to a ket in a well defined way.");
    }
    return kets[ket_of_states[index_state]];
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::begin() const {
    return kets.begin();
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::end() const {
    return kets.end();
}

template <typename Derived>
Basis<Derived>::Iterator::Iterator(typename ketvec_t::const_iterator it) : it{std::move(it)} {}

template <typename Derived>
bool Basis<Derived>::Iterator::operator!=(const Iterator &other) const {
    return other.it != it;
}

template <typename Derived>
std::shared_ptr<const typename Basis<Derived>::ket_t> Basis<Derived>::Iterator::operator*() const {
    return *it;
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++it;
    return *this;
}

template <typename Derived>
size_t Basis<Derived>::get_number_of_states() const {
    return coefficients.matrix.cols();
}

template <typename Derived>
size_t Basis<Derived>::get_number_of_kets() const {
    return coefficients.matrix.rows();
}

template <typename Derived>
const Transformation<typename Basis<Derived>::scalar_t> &
Basis<Derived>::get_transformation() const {
    return coefficients;
}

template <typename Derived>
Transformation<typename Basis<Derived>::scalar_t>
Basis<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    Transformation<scalar_t> transformation{{static_cast<Eigen::Index>(coefficients.matrix.rows()),
                                             static_cast<Eigen::Index>(coefficients.matrix.rows())},
                                            {TransformationType::ROTATE}};

    std::vector<Eigen::Triplet<scalar_t>> entries;

    for (size_t idx_initial = 0; idx_initial < kets.size(); ++idx_initial) {
        real_t f = kets[idx_initial]->get_quantum_number_f();
        real_t m_initial = kets[idx_initial]->get_quantum_number_m();
        for (real_t m_final = -f; m_final <= f; ++m_final) {
            auto val = wigner::wigner_uppercase_d_matrix<scalar_t>(f, m_initial, m_final, alpha,
                                                                   beta, gamma);
            size_t idx_final = ket_id_to_index.at(
                kets[idx_initial]->get_id_for_different_quantum_number_m(m_final));
            entries.emplace_back(idx_final, idx_initial, val);
        }
    }

    transformation.matrix.setFromTriplets(entries.begin(), entries.end());
    transformation.matrix.makeCompressed();

    return transformation;
}

template <typename Derived>
Sorting Basis<Derived>::get_sorter(const std::vector<TransformationType> &labels) const {
    perform_sorter_checks(labels);

    // Throw a meaningful error if sorting by energy is requested as this might be a common mistake
    if (std::find(labels.begin(), labels.end(), TransformationType::SORT_BY_ENERGY) !=
        labels.end()) {
        throw std::invalid_argument("States do not store the energy and thus can not be sorted by "
                                    "the energy. Use an energy operator instead.");
    }

    // Initialize transformation
    Sorting transformation;
    transformation.matrix.resize(coefficients.matrix.cols());
    transformation.matrix.setIdentity();

    // Get the sorter
    get_sorter_without_checks(labels, transformation);

    // Check if all labels have been used for sorting
    if (labels != transformation.transformation_type) {
        throw std::invalid_argument("The states could not be sorted by all the requested labels.");
    }

    return transformation;
}

template <typename Derived>
std::vector<IndicesOfBlock>
Basis<Derived>::get_indices_of_blocks(const std::vector<TransformationType> &labels) const {
    perform_sorter_checks(labels);

    std::set<TransformationType> unique_labels(labels.begin(), labels.end());
    perform_blocks_checks(unique_labels);

    // Get the blocks
    IndicesOfBlocksCreator blocks_creator({0, static_cast<size_t>(coefficients.matrix.cols())});
    get_indices_of_blocks_without_checks(unique_labels, blocks_creator);

    return blocks_creator.create();
}

template <typename Derived>
void Basis<Derived>::get_sorter_without_checks(const std::vector<TransformationType> &labels,
                                               Sorting &transformation) const {
    int *perm_begin = transformation.matrix.indices().data();
    int *perm_end = perm_begin + coefficients.matrix.cols();
    int *perm_back = perm_end - 1;

    // Sort the vector based on the requested labels
    std::stable_sort(perm_begin, perm_end, [&](int a, int b) {
        for (const auto &label : labels) {
            switch (label) {
            case TransformationType::SORT_BY_PARITY:
                if (parity_of_states[a] != parity_of_states[b]) {
                    return parity_of_states[a] < parity_of_states[b];
                }
                break;
            case TransformationType::SORT_BY_QUANTUM_NUMBER_M:
                if (quantum_number_m_of_states[a] != quantum_number_m_of_states[b]) {
                    return quantum_number_m_of_states[a] < quantum_number_m_of_states[b];
                }
                break;
            case TransformationType::SORT_BY_QUANTUM_NUMBER_F:
                if (quantum_number_f_of_states[a] != quantum_number_f_of_states[b]) {
                    return quantum_number_f_of_states[a] < quantum_number_f_of_states[b];
                }
                break;
            case TransformationType::SORT_BY_KET:
                if (ket_of_states[a] != ket_of_states[b]) {
                    return ket_of_states[a] < ket_of_states[b];
                }
                break;
            default:
                std::abort(); // Can't happen because of previous checks
            }
        }
        return false; // Elements are equal
    });

    // Check for invalid values and add transformation types
    for (const auto &label : labels) {
        switch (label) {
        case TransformationType::SORT_BY_PARITY:
            if (parity_of_states[*perm_back] == Parity::UNKNOWN) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by the parity.");
            }
            transformation.transformation_type.push_back(TransformationType::SORT_BY_PARITY);
            break;
        case TransformationType::SORT_BY_QUANTUM_NUMBER_M:
            if (quantum_number_m_of_states[*perm_back] == std::numeric_limits<real_t>::max()) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by the quantum number m.");
            }
            transformation.transformation_type.push_back(
                TransformationType::SORT_BY_QUANTUM_NUMBER_M);
            break;
        case TransformationType::SORT_BY_QUANTUM_NUMBER_F:
            if (quantum_number_f_of_states[*perm_back] == std::numeric_limits<real_t>::max()) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by the quantum number f.");
            }
            transformation.transformation_type.push_back(
                TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            break;
        case TransformationType::SORT_BY_KET:
            if (ket_of_states[*perm_back] == std::numeric_limits<int>::max()) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by kets.");
            }
            transformation.transformation_type.push_back(TransformationType::SORT_BY_KET);
            break;
        default:
            std::abort(); // Can't happen because of previous checks
        }
    }
}

template <typename Derived>
void Basis<Derived>::get_indices_of_blocks_without_checks(
    const std::set<TransformationType> &unique_labels,
    IndicesOfBlocksCreator &blocks_creator) const {
    auto last_quantum_number_f = quantum_number_f_of_states[0];
    auto last_quantum_number_m = quantum_number_m_of_states[0];
    auto last_parity = parity_of_states[0];
    auto last_ket = ket_of_states[0];

    for (int i = 0; i < coefficients.matrix.cols(); ++i) {
        for (auto label : unique_labels) {
            if (label == TransformationType::SORT_BY_QUANTUM_NUMBER_F) {
                if (quantum_number_f_of_states[i] != last_quantum_number_f) {
                    blocks_creator.add(i);
                    break;
                }
            } else if (label == TransformationType::SORT_BY_QUANTUM_NUMBER_M) {
                if (quantum_number_m_of_states[i] != last_quantum_number_m) {
                    blocks_creator.add(i);
                    break;
                }
            } else if (label == TransformationType::SORT_BY_PARITY) {
                if (parity_of_states[i] != last_parity) {
                    blocks_creator.add(i);
                    break;
                }
            } else if (label == TransformationType::SORT_BY_KET) {
                if (ket_of_states[i] != last_ket) {
                    blocks_creator.add(i);
                    break;
                }
            }
        }
        last_quantum_number_f = quantum_number_f_of_states[i];
        last_quantum_number_m = quantum_number_m_of_states[i];
        last_parity = parity_of_states[i];
        last_ket = ket_of_states[i];
    }
}

template <typename Derived>
std::shared_ptr<Derived> Basis<Derived>::transformed(const Sorting &transformation) const {
    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    // Apply the transformation
    transformed->coefficients.matrix = transformed->coefficients.matrix * transformation.matrix;
    transformed->coefficients.transformation_type = transformation.transformation_type;

    transformed->quantum_number_f_of_states.resize(transformation.matrix.size());
    transformed->quantum_number_m_of_states.resize(transformation.matrix.size());
    transformed->parity_of_states.resize(transformation.matrix.size());
    transformed->ket_of_states.resize(transformation.matrix.size());

    for (int i = 0; i < transformation.matrix.size(); ++i) {
        transformed->quantum_number_f_of_states[i] =
            quantum_number_f_of_states[transformation.matrix.indices()[i]];
        transformed->quantum_number_m_of_states[i] =
            quantum_number_m_of_states[transformation.matrix.indices()[i]];
        transformed->parity_of_states[i] = parity_of_states[transformation.matrix.indices()[i]];
        transformed->ket_of_states[i] = ket_of_states[transformation.matrix.indices()[i]];
    }

    return transformed;
}

template <typename Derived>
std::shared_ptr<Derived>
Basis<Derived>::transformed(const Transformation<scalar_t> &transformation) const {
    // If the transformation is a rotation, it should be a rotation and nothing else
    bool is_rotation = false;
    for (auto t : transformation.transformation_type) {
        if (t == TransformationType::ROTATE) {
            is_rotation = true;
            break;
        }
    }
    if (is_rotation && transformation.transformation_type.size() != 1) {
        throw std::invalid_argument("A rotation can not be combined with other transformations.");
    }

    // To apply a rotation, the object must only be sorted but other transformations are not allowed
    if (is_rotation) {
        for (auto t : coefficients.transformation_type) {
            if (!utils::is_sorting(t)) {
                throw std::runtime_error(
                    "If the object was transformed by a different transformation "
                    "than sorting, it can not be rotated.");
            }
        }
    }

    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    // Apply the transformation
    transformed->coefficients.matrix = coefficients.matrix * transformation.matrix;
    transformed->coefficients.transformation_type = transformation.transformation_type;

    Eigen::SparseMatrix<real_t> probs = transformation.matrix.cwiseAbs2().transpose();

    {
        auto map = Eigen::Map<const Eigen::VectorX<real_t>>(quantum_number_f_of_states.data(),
                                                            quantum_number_f_of_states.size());
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        Eigen::VectorX<real_t> diff = (val.cwiseAbs2() - sq).cwiseAbs();
        transformed->quantum_number_f_of_states.resize(probs.rows());

        for (size_t i = 0; i < transformed->quantum_number_f_of_states.size(); ++i) {
            if (diff[i] < 10 * std::numeric_limits<real_t>::epsilon()) {
                transformed->quantum_number_f_of_states[i] = val[i];
            } else {
                transformed->quantum_number_f_of_states[i] = std::numeric_limits<real_t>::max();
            }
        }
    }

    {
        auto map = Eigen::Map<const Eigen::VectorX<real_t>>(quantum_number_m_of_states.data(),
                                                            quantum_number_m_of_states.size());
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        Eigen::VectorX<real_t> diff = (val.cwiseAbs2() - sq).cwiseAbs();
        transformed->quantum_number_m_of_states.resize(probs.rows());

        for (size_t i = 0; i < transformed->quantum_number_m_of_states.size(); ++i) {
            if (diff[i] < 10 * std::numeric_limits<real_t>::epsilon()) {
                transformed->quantum_number_m_of_states[i] = val[i];
            } else {
                transformed->quantum_number_m_of_states[i] = std::numeric_limits<real_t>::max();
            }
        }
    }

    {
        using utype = std::underlying_type<Parity>::type;
        Eigen::VectorX<real_t> map(parity_of_states.size());
        for (size_t i = 0; i < parity_of_states.size(); ++i) {
            map[i] = static_cast<utype>(parity_of_states[i]);
        }
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        Eigen::VectorX<real_t> diff = (val.cwiseAbs2() - sq).cwiseAbs();
        transformed->parity_of_states.resize(probs.rows());

        for (size_t i = 0; i < transformed->parity_of_states.size(); ++i) {
            if (diff[i] < 10 * std::numeric_limits<real_t>::epsilon()) {
                transformed->parity_of_states[i] = static_cast<Parity>(std::lround(val[i]));
            } else {
                transformed->parity_of_states[i] = Parity::UNKNOWN;
            }
        }
    }

    {
        transformed->ket_of_states.resize(transformed->coefficients.matrix.cols());
        std::fill(transformed->ket_of_states.begin(), transformed->ket_of_states.end(),
                  std::numeric_limits<int>::max());
        std::vector<real_t> map_idx_to_max(transformed->coefficients.matrix.cols(), 0);
        for (int row = 0; row < transformed->coefficients.matrix.outerSize(); ++row) {
            for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(
                     transformed->coefficients.matrix, row);
                 it; ++it) {
                if (std::abs(it.value()) > map_idx_to_max[it.col()]) {
                    map_idx_to_max[it.col()] = std::abs(it.value());
                    transformed->ket_of_states[it.col()] = row;
                }
            }
        }
    }

    return transformed;
}

// Explicit instantiations
template class Basis<BasisAtom<float>>;
template class Basis<BasisAtom<double>>;
template class Basis<BasisAtom<std::complex<float>>>;
template class Basis<BasisAtom<std::complex<double>>>;
template class Basis<BasisClassicalLight<float>>;
template class Basis<BasisClassicalLight<double>>;
template class Basis<BasisClassicalLight<std::complex<float>>>;
template class Basis<BasisClassicalLight<std::complex<double>>>;
} // namespace pairinteraction
