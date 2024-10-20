#include "pairinteraction/basis/Basis.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisClassicalLight.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetClassicalLight.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/wigner.hpp"

#include <numeric>
#include <set>

namespace pairinteraction {

template <typename Scalar>
class BasisAtom;
template <typename Scalar>
class BasisClassicalLight;

template <typename Derived>
const std::string &Basis<Derived>::get_id_of_kets() const {
    return id_of_kets;
}

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
Basis<Derived>::Basis(ketvec_t &&kets, std::string &&id_of_kets)
    : kets(std::move(kets)), coefficients{{static_cast<Eigen::Index>(this->kets.size()),
                                           static_cast<Eigen::Index>(this->kets.size())},
                                          {TransformationType::SORT_BY_KET}},
      id_of_kets(std::move(id_of_kets)) {
    if (this->kets.empty()) {
        throw std::invalid_argument("The basis must contain at least one element.");
    }
    state_index_to_quantum_number_f.reserve(this->kets.size());
    state_index_to_quantum_number_m.reserve(this->kets.size());
    state_index_to_parity.reserve(this->kets.size());
    ket_id_to_ket_index.reserve(this->kets.size());
    size_t index = 0;
    for (const auto &ket : this->kets) {
        state_index_to_quantum_number_f.push_back(ket->get_quantum_number_f());
        state_index_to_quantum_number_m.push_back(ket->get_quantum_number_m());
        state_index_to_parity.push_back(ket->get_parity());
        ket_id_to_ket_index[ket->get_id()] = index++;
        if (ket->get_quantum_number_f() == std::numeric_limits<real_t>::max()) {
            _has_quantum_number_f = false;
        }
        if (ket->get_quantum_number_m() == std::numeric_limits<real_t>::max()) {
            _has_quantum_number_m = false;
        }
        if (ket->get_parity() == Parity::UNKNOWN) {
            _has_parity = false;
        }
    }
    state_index_to_ket_index.resize(this->kets.size());
    std::iota(state_index_to_ket_index.begin(), state_index_to_ket_index.end(), 0);
    ket_index_to_state_index.resize(this->kets.size());
    std::iota(ket_index_to_state_index.begin(), ket_index_to_state_index.end(), 0);
    coefficients.matrix.setIdentity();
}

template <typename Derived>
bool Basis<Derived>::has_quantum_number_f() const {
    return _has_quantum_number_f;
}

template <typename Derived>
bool Basis<Derived>::has_quantum_number_m() const {
    return _has_quantum_number_m;
}

template <typename Derived>
bool Basis<Derived>::has_parity() const {
    return _has_parity;
}

template <typename Derived>
bool Basis<Derived>::has_ket_index() const {
    return _has_ket_index;
}

template <typename Derived>
bool Basis<Derived>::has_ket_index(size_t ket_id) const {
    return ket_id_to_ket_index.count(ket_id) > 0;
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
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t, Eigen::RowMajor> &
Basis<Derived>::get_coefficients() {
    return coefficients.matrix;
}

template <typename Derived>
size_t Basis<Derived>::get_ket_index(size_t ket_id) const {
    return ket_id_to_ket_index.at(ket_id);
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t, Eigen::RowMajor>
Basis<Derived>::get_amplitudes(std::shared_ptr<const ket_t> ket) const {
    if (!has_ket_index(ket->get_id())) {
        throw std::invalid_argument("The ket does not belong to the basis.");
    }
    // The following line is a more efficient alternative to
    // "get_amplitudes(get_state_from_ket(ket))"
    return coefficients.matrix.row(ket_id_to_ket_index.at(ket->get_id()));
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t, Eigen::RowMajor>
Basis<Derived>::get_amplitudes(std::shared_ptr<const Derived> other) const {
    if (other->get_id_of_kets() != get_id_of_kets()) {
        throw std::invalid_argument("The other object must be expressed using the same kets.");
    }
    return other->coefficients.matrix.adjoint() * coefficients.matrix;
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::real_t, Eigen::RowMajor>
Basis<Derived>::get_overlaps(std::shared_ptr<const ket_t> ket) const {
    return get_amplitudes(ket).cwiseAbs2();
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::real_t, Eigen::RowMajor>
Basis<Derived>::get_overlaps(std::shared_ptr<const Derived> other) const {
    return get_amplitudes(other).cwiseAbs2();
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_quantum_number_f(size_t state_index) const {
    real_t quantum_number_f = state_index_to_quantum_number_f.at(state_index);
    if (quantum_number_f == std::numeric_limits<real_t>::max()) {
        throw std::invalid_argument("The state does not have a well-defined quantum number f.");
    }
    return quantum_number_f;
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_quantum_number_m(size_t state_index) const {
    real_t quantum_number_m = state_index_to_quantum_number_m.at(state_index);
    if (quantum_number_m == std::numeric_limits<real_t>::max()) {
        throw std::invalid_argument("The state does not have a well-defined quantum number m.");
    }
    return quantum_number_m;
}

template <typename Derived>
Parity Basis<Derived>::get_parity(size_t state_index) const {
    Parity parity = state_index_to_parity.at(state_index);
    if (parity == Parity::UNKNOWN) {
        throw std::invalid_argument("The state does not have a well-defined parity.");
    }
    return parity;
}

template <typename Derived>
std::shared_ptr<const typename Basis<Derived>::ket_t>
Basis<Derived>::get_ket_with_largest_overlap(size_t state_index) const {
    size_t ket_index = state_index_to_ket_index.at(state_index);
    if (ket_index == std::numeric_limits<int>::max()) {
        throw std::invalid_argument("The state does not belong to a ket in a well-defined way.");
    }
    return kets[ket_index];
}

template <typename Derived>
std::shared_ptr<const Derived>
Basis<Derived>::get_state_with_largest_overlap(size_t ket_index) const {
    int state_index = ket_index_to_state_index[ket_index];
    if (state_index == std::numeric_limits<int>::max()) {
        throw std::runtime_error("The ket does not belong to a state in a well-defined way.");
    }
    // Create a copy of the current object
    auto restricted = std::make_shared<Derived>(derived());

    // Restrict the copy to the state with the largest overlap
    restricted->coefficients.matrix = restricted->coefficients.matrix.col(state_index);

    std::fill(restricted->ket_index_to_state_index.begin(),
              restricted->ket_index_to_state_index.end(), std::numeric_limits<int>::max());
    restricted->ket_index_to_state_index[ket_index] = 0;

    restricted->state_index_to_quantum_number_f = {state_index_to_quantum_number_f[state_index]};
    restricted->state_index_to_quantum_number_m = {state_index_to_quantum_number_m[state_index]};
    restricted->state_index_to_parity = {state_index_to_parity[state_index]};
    restricted->state_index_to_ket_index = {state_index_to_ket_index[state_index]};

    restricted->_has_quantum_number_f =
        restricted->state_index_to_quantum_number_f[0] != std::numeric_limits<real_t>::max();
    restricted->_has_quantum_number_m =
        restricted->state_index_to_quantum_number_m[0] != std::numeric_limits<real_t>::max();
    restricted->_has_parity = restricted->state_index_to_parity[0] != Parity::UNKNOWN;
    restricted->_has_ket_index =
        restricted->state_index_to_ket_index[0] != std::numeric_limits<int>::max();

    return restricted;
}

template <typename Derived>
std::shared_ptr<const Derived>
Basis<Derived>::get_state_with_largest_overlap(std::shared_ptr<const ket_t> ket) const {
    if (!has_ket_index(ket->get_id())) {
        throw std::invalid_argument("The ket does not belong to the basis.");
    }
    return get_state_with_largest_overlap(ket_id_to_ket_index.at(ket->get_id()));
}

template <typename Derived>
std::shared_ptr<const Derived> Basis<Derived>::get_state_from_ket(size_t ket_index) const {
    // Create a copy of the current object
    auto created = std::make_shared<Derived>(derived());

    // Fill the copy with the state corresponding to the ket index
    created->coefficients.matrix =
        Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>(coefficients.matrix.rows(), 1);
    created->coefficients.matrix.coeffRef(ket_index, 0) = 1;
    created->coefficients.matrix.makeCompressed();

    std::fill(created->ket_index_to_state_index.begin(), created->ket_index_to_state_index.end(),
              std::numeric_limits<int>::max());
    created->ket_index_to_state_index[ket_index] = 0;

    created->state_index_to_quantum_number_f = {kets[ket_index]->get_quantum_number_f()};
    created->state_index_to_quantum_number_m = {kets[ket_index]->get_quantum_number_m()};
    created->state_index_to_parity = {kets[ket_index]->get_parity()};
    created->state_index_to_ket_index = {ket_index};

    created->_has_quantum_number_f =
        created->state_index_to_quantum_number_f[0] != std::numeric_limits<real_t>::max();
    created->_has_quantum_number_m =
        created->state_index_to_quantum_number_m[0] != std::numeric_limits<real_t>::max();
    created->_has_parity = created->state_index_to_parity[0] != Parity::UNKNOWN;
    created->_has_ket_index =
        created->state_index_to_ket_index[0] != std::numeric_limits<int>::max();

    return created;
}

template <typename Derived>
std::shared_ptr<const Derived>
Basis<Derived>::get_state_from_ket(std::shared_ptr<const ket_t> ket) const {
    if (!has_ket_index(ket->get_id())) {
        throw std::invalid_argument("The ket does not belong to the basis.");
    }
    return get_state_from_ket(ket_id_to_ket_index.at(ket->get_id()));
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
            size_t idx_final = ket_id_to_ket_index.at(
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
    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    int *perm_begin = transformation.matrix.indices().data();
    int *perm_end = perm_begin + coefficients.matrix.cols();
    int *perm_back = perm_end - 1;

    // Sort the vector based on the requested labels
    std::stable_sort(perm_begin, perm_end, [&](int a, int b) {
        for (const auto &label : labels) {
            switch (label) {
            case TransformationType::SORT_BY_PARITY:
                if (state_index_to_parity[a] != state_index_to_parity[b]) {
                    return state_index_to_parity[a] < state_index_to_parity[b];
                }
                break;
            case TransformationType::SORT_BY_QUANTUM_NUMBER_M:
                if (std::abs(state_index_to_quantum_number_m[a] -
                             state_index_to_quantum_number_m[b]) > numerical_precision) {
                    return state_index_to_quantum_number_m[a] < state_index_to_quantum_number_m[b];
                }
                break;
            case TransformationType::SORT_BY_QUANTUM_NUMBER_F:
                if (std::abs(state_index_to_quantum_number_f[a] -
                             state_index_to_quantum_number_f[b]) > numerical_precision) {
                    return state_index_to_quantum_number_f[a] < state_index_to_quantum_number_f[b];
                }
                break;
            case TransformationType::SORT_BY_KET:
                if (state_index_to_ket_index[a] != state_index_to_ket_index[b]) {
                    return state_index_to_ket_index[a] < state_index_to_ket_index[b];
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
            if (state_index_to_parity[*perm_back] == Parity::UNKNOWN) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by the parity.");
            }
            transformation.transformation_type.push_back(TransformationType::SORT_BY_PARITY);
            break;
        case TransformationType::SORT_BY_QUANTUM_NUMBER_M:
            if (state_index_to_quantum_number_m[*perm_back] == std::numeric_limits<real_t>::max()) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by the quantum number m.");
            }
            transformation.transformation_type.push_back(
                TransformationType::SORT_BY_QUANTUM_NUMBER_M);
            break;
        case TransformationType::SORT_BY_QUANTUM_NUMBER_F:
            if (state_index_to_quantum_number_f[*perm_back] == std::numeric_limits<real_t>::max()) {
                throw std::invalid_argument(
                    "States cannot be labeled and thus not sorted by the quantum number f.");
            }
            transformation.transformation_type.push_back(
                TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            break;
        case TransformationType::SORT_BY_KET:
            if (state_index_to_ket_index[*perm_back] == std::numeric_limits<int>::max()) {
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
    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    auto last_quantum_number_f = state_index_to_quantum_number_f[0];
    auto last_quantum_number_m = state_index_to_quantum_number_m[0];
    auto last_parity = state_index_to_parity[0];
    auto last_ket = state_index_to_ket_index[0];

    for (int i = 0; i < coefficients.matrix.cols(); ++i) {
        for (auto label : unique_labels) {
            if (label == TransformationType::SORT_BY_QUANTUM_NUMBER_F) {
                if (std::abs(state_index_to_quantum_number_f[i] - last_quantum_number_f) >
                    numerical_precision) {
                    blocks_creator.add(i);
                    break;
                }
            } else if (label == TransformationType::SORT_BY_QUANTUM_NUMBER_M) {
                if (std::abs(state_index_to_quantum_number_m[i] - last_quantum_number_m) >
                    numerical_precision) {
                    blocks_creator.add(i);
                    break;
                }
            } else if (label == TransformationType::SORT_BY_PARITY) {
                if (state_index_to_parity[i] != last_parity) {
                    blocks_creator.add(i);
                    break;
                }
            } else if (label == TransformationType::SORT_BY_KET) {
                if (state_index_to_ket_index[i] != last_ket) {
                    blocks_creator.add(i);
                    break;
                }
            }
        }
        last_quantum_number_f = state_index_to_quantum_number_f[i];
        last_quantum_number_m = state_index_to_quantum_number_m[i];
        last_parity = state_index_to_parity[i];
        last_ket = state_index_to_ket_index[i];
    }
}

template <typename Derived>
std::shared_ptr<const Derived> Basis<Derived>::transformed(const Sorting &transformation) const {
    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    // Apply the transformation
    transformed->coefficients.matrix = transformed->coefficients.matrix * transformation.matrix;
    transformed->coefficients.transformation_type = transformation.transformation_type;

    transformed->state_index_to_quantum_number_f.resize(transformation.matrix.size());
    transformed->state_index_to_quantum_number_m.resize(transformation.matrix.size());
    transformed->state_index_to_parity.resize(transformation.matrix.size());
    transformed->state_index_to_ket_index.resize(transformation.matrix.size());
    transformed->ket_index_to_state_index.resize(transformation.matrix.size());

    for (int i = 0; i < transformation.matrix.size(); ++i) {
        transformed->state_index_to_quantum_number_f[i] =
            state_index_to_quantum_number_f[transformation.matrix.indices()[i]];
        transformed->state_index_to_quantum_number_m[i] =
            state_index_to_quantum_number_m[transformation.matrix.indices()[i]];
        transformed->state_index_to_parity[i] =
            state_index_to_parity[transformation.matrix.indices()[i]];
        transformed->state_index_to_ket_index[i] =
            state_index_to_ket_index[transformation.matrix.indices()[i]];
        transformed->ket_index_to_state_index[transformation.matrix.indices()[i]] = i;
    }

    return transformed;
}

template <typename Derived>
std::shared_ptr<const Derived>
Basis<Derived>::transformed(const Transformation<scalar_t> &transformation) const {
    real_t numerical_precision =
        10 * std::sqrt(coefficients.matrix.rows()) * std::numeric_limits<real_t>::epsilon();

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
    // If a quantum number turns out to be conserved by the transformation, it will be
    // rounded to the nearest half integer to avoid loss of numerical_precision.
    transformed->coefficients.matrix = coefficients.matrix * transformation.matrix;
    transformed->coefficients.transformation_type = transformation.transformation_type;

    Eigen::SparseMatrix<real_t> probs = transformation.matrix.cwiseAbs2().transpose();

    {
        auto map = Eigen::Map<const Eigen::VectorX<real_t>>(state_index_to_quantum_number_f.data(),
                                                            state_index_to_quantum_number_f.size());
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        Eigen::VectorX<real_t> diff = (val.cwiseAbs2() - sq).cwiseAbs();
        transformed->state_index_to_quantum_number_f.resize(probs.rows());

        for (size_t i = 0; i < transformed->state_index_to_quantum_number_f.size(); ++i) {
            if (diff[i] < numerical_precision) {
                transformed->state_index_to_quantum_number_f[i] = std::round(val[i] * 2) / 2;
            } else {
                transformed->state_index_to_quantum_number_f[i] =
                    std::numeric_limits<real_t>::max();
                transformed->_has_quantum_number_f = false;
            }
        }
    }

    {
        auto map = Eigen::Map<const Eigen::VectorX<real_t>>(state_index_to_quantum_number_m.data(),
                                                            state_index_to_quantum_number_m.size());
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        Eigen::VectorX<real_t> diff = (val.cwiseAbs2() - sq).cwiseAbs();
        transformed->state_index_to_quantum_number_m.resize(probs.rows());

        for (size_t i = 0; i < transformed->state_index_to_quantum_number_m.size(); ++i) {
            if (diff[i] < numerical_precision) {
                transformed->state_index_to_quantum_number_m[i] = std::round(val[i] * 2) / 2;
            } else {
                transformed->state_index_to_quantum_number_m[i] =
                    std::numeric_limits<real_t>::max();
                transformed->_has_quantum_number_m = false;
            }
        }
    }

    {
        using utype = std::underlying_type<Parity>::type;
        Eigen::VectorX<real_t> map(state_index_to_parity.size());
        for (size_t i = 0; i < state_index_to_parity.size(); ++i) {
            map[i] = static_cast<utype>(state_index_to_parity[i]);
        }
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        Eigen::VectorX<real_t> diff = (val.cwiseAbs2() - sq).cwiseAbs();
        transformed->state_index_to_parity.resize(probs.rows());

        for (size_t i = 0; i < transformed->state_index_to_parity.size(); ++i) {
            if (diff[i] < numerical_precision) {
                transformed->state_index_to_parity[i] = static_cast<Parity>(std::lround(val[i]));
            } else {
                transformed->state_index_to_parity[i] = Parity::UNKNOWN;
                transformed->_has_parity = false;
            }
        }
    }

    {
        transformed->state_index_to_ket_index.resize(transformed->coefficients.matrix.cols());
        transformed->ket_index_to_state_index.resize(transformed->coefficients.matrix.rows());
        std::fill(transformed->state_index_to_ket_index.begin(),
                  transformed->state_index_to_ket_index.end(), std::numeric_limits<int>::max());
        std::fill(transformed->ket_index_to_state_index.begin(),
                  transformed->ket_index_to_state_index.end(), std::numeric_limits<int>::max());
        std::vector<real_t> map_idx_to_max(transformed->coefficients.matrix.cols(), 0);
        for (int row = 0; row < transformed->coefficients.matrix.outerSize(); ++row) {
            for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(
                     transformed->coefficients.matrix, row);
                 it; ++it) {
                if (std::abs(it.value()) > map_idx_to_max[it.col()]) {
                    map_idx_to_max[it.col()] = std::abs(it.value());
                    transformed->state_index_to_ket_index[it.col()] = row;
                    transformed->ket_index_to_state_index[row] = it.col();
                }
            }
        }

        // Set _has_ket_index to false if any of the ket indices is invalid
        transformed->_has_ket_index = true;
        for (auto ket_index : transformed->state_index_to_ket_index) {
            if (ket_index == std::numeric_limits<int>::max()) {
                transformed->_has_ket_index = false;
                break;
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
template class Basis<BasisCombined<float>>;
template class Basis<BasisCombined<double>>;
template class Basis<BasisCombined<std::complex<float>>>;
template class Basis<BasisCombined<std::complex<double>>>;
} // namespace pairinteraction
