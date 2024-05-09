#include "basis/Basis.hpp"
#include "ket/Ket.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetClassicalLight.hpp"
#include "utils/euler.hpp"
#include "utils/wigner.hpp"

#include <numeric>
#include <set>

template <typename Derived>
Basis<Derived>::Basis(ketvec_t &&kets)
    : kets(std::move(kets)), coefficients{{static_cast<Eigen::Index>(this->kets.size()),
                                           static_cast<Eigen::Index>(this->kets.size())},
                                          {TransformationType::SORT_BY_KET}} {
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
float Basis<Derived>::get_quantum_number_f(size_t index_state) const {
    return quantum_number_f_of_states[index_state];
}

template <typename Derived>
float Basis<Derived>::get_quantum_number_m(size_t index_state) const {
    return quantum_number_m_of_states[index_state];
}

template <typename Derived>
int Basis<Derived>::get_parity(size_t index_state) const {
    return parity_of_states[index_state];
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::begin() const {
    return Iterator(*this, 0);
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::end() const {
    return Iterator(*this, kets.size());
}

template <typename Derived>
Basis<Derived>::Iterator::Iterator(const Basis<Derived> &basis, size_t index)
    : basis(basis), index(index) {}

template <typename Derived>
bool Basis<Derived>::Iterator::operator!=(const Iterator &other) const {
    return index != other.index;
}

template <typename Derived>
std::shared_ptr<const typename Basis<Derived>::ket_t> Basis<Derived>::Iterator::operator*() const {
    return basis.kets[index];
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++index;
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
        float f = kets[idx_initial]->get_quantum_number_f();
        float m_initial = kets[idx_initial]->get_quantum_number_m();
        for (float m_final = -f; m_final <= f; ++m_final) {
            scalar_t val = wigner::wigner_uppercase_d_matrix<scalar_t>(f, m_initial, m_final, alpha,
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
Sorting Basis<Derived>::get_sorter(TransformationType label) const {
    // Check that the label is a valid sorting label
    if (!utils::is_sorting(label)) {
        throw std::invalid_argument("The label is not a valid sorting label.");
    }

    // Get the sorter
    Sorting transformation = get_sorter_without_checks(label);

    // Throw a meaningful error if sorting by energy is requested as this might be a common mistake
    if ((label & TransformationType::SORT_BY_ENERGY) == TransformationType::SORT_BY_ENERGY) {
        throw std::invalid_argument("States do not store the energy and thus can not be sorted by "
                                    "the energy. Use an energy operator instead.");
    }

    // Check if the full label has been used for sorting
    if (!utils::is_comprised_by_label(label, transformation.transformation_type)) {
        throw std::invalid_argument("The states could not be sorted by the requested label.");
    }

    return transformation;
}

template <typename Derived>
Blocks Basis<Derived>::get_blocks(TransformationType label) const {
    // Check that the label is a valid sorting label
    if (!utils::is_sorting(label)) {
        throw std::invalid_argument("The label is not a valid sorting label.");
    }

    // Check if the states are sorted by the requested label
    if (!utils::is_sorted_by_label(label, get_transformation().transformation_type)) {
        throw std::invalid_argument("The states are not sorted by the requested label.");
    }

    // Get the blocks
    Blocks blocks = get_blocks_without_checks(label);

    // Throw a meaningful error if getting the blocks by energy is requested as this might be a
    // common mistake
    if ((label & TransformationType::SORT_BY_ENERGY) == TransformationType::SORT_BY_ENERGY) {
        throw std::invalid_argument("States do not store the energy and thus no energy blocks can "
                                    "be obtained. Use an energy operator instead.");
    }

    // Check if the full label has been used for getting the blocks
    if (!utils::is_comprised_by_label(label, blocks.transformation_type)) {
        throw std::invalid_argument("The blocks could not be obtained by the requested label.");
    }

    return blocks;
}

template <typename Derived>
Sorting Basis<Derived>::get_sorter_without_checks(TransformationType label) const {
    Sorting transformation;

    std::vector<int> perm(coefficients.matrix.cols());
    std::iota(perm.begin(), perm.end(), 0);

    if ((label & TransformationType::SORT_BY_QUANTUM_NUMBER_F) ==
        TransformationType::SORT_BY_QUANTUM_NUMBER_F) {
        std::stable_sort(perm.begin(), perm.end(), [&](int i, int j) {
            return quantum_number_f_of_states[i] < quantum_number_f_of_states[j];
        });

        if (quantum_number_f_of_states[perm.back()] == std::numeric_limits<float>::max()) {
            throw std::invalid_argument(
                "States cannot be labeled and thus not sorted by the quantum number f.");
        }

        transformation.transformation_type.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
    }

    if ((label & TransformationType::SORT_BY_QUANTUM_NUMBER_M) ==
        TransformationType::SORT_BY_QUANTUM_NUMBER_M) {
        std::stable_sort(perm.begin(), perm.end(), [&](int i, int j) {
            return quantum_number_m_of_states[i] < quantum_number_m_of_states[j];
        });

        if (quantum_number_m_of_states[perm.back()] == std::numeric_limits<float>::max()) {
            throw std::invalid_argument(
                "States cannot be labeled and thus not sorted by the quantum number m.");
        }

        transformation.transformation_type.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }

    if ((label & TransformationType::SORT_BY_PARITY) == TransformationType::SORT_BY_PARITY) {
        std::stable_sort(perm.begin(), perm.end(),
                         [&](int i, int j) { return parity_of_states[i] < parity_of_states[j]; });

        if (parity_of_states[perm.back()] == std::numeric_limits<int>::max()) {
            throw std::invalid_argument(
                "States cannot be labeled and thus not sorted by the parity.");
        }

        transformation.transformation_type.push_back(TransformationType::SORT_BY_PARITY);
    }

    if ((label & TransformationType::SORT_BY_KET) == TransformationType::SORT_BY_KET) {
        std::stable_sort(perm.begin(), perm.end(),
                         [&](int i, int j) { return ket_of_states[i] < ket_of_states[j]; });

        transformation.transformation_type.push_back(TransformationType::SORT_BY_KET);
    }

    transformation.matrix.indices() = Eigen::Map<const Eigen::VectorXi>(perm.data(), perm.size());

    return transformation;
}

template <typename Derived>
Blocks Basis<Derived>::get_blocks_without_checks(TransformationType label) const {
    Blocks blocks;

    float last_quantum_number_f = quantum_number_f_of_states[0];
    float last_quantum_number_m = quantum_number_m_of_states[0];
    int last_parity = parity_of_states[0];
    int last_ket = ket_of_states[0];
    blocks.start.push_back(0);

    for (int i = 0; i < coefficients.matrix.cols(); ++i) {
        if ((label & TransformationType::SORT_BY_QUANTUM_NUMBER_F) ==
                TransformationType::SORT_BY_QUANTUM_NUMBER_F &&
            quantum_number_f_of_states[i] != last_quantum_number_f) {
            blocks.start.push_back(i);
        } else if ((label & TransformationType::SORT_BY_QUANTUM_NUMBER_M) ==
                       TransformationType::SORT_BY_QUANTUM_NUMBER_M &&
                   quantum_number_m_of_states[i] != last_quantum_number_m) {
            blocks.start.push_back(i);
        } else if ((label & TransformationType::SORT_BY_PARITY) ==
                       TransformationType::SORT_BY_PARITY &&
                   parity_of_states[i] != last_parity) {
            blocks.start.push_back(i);
        } else if ((label & TransformationType::SORT_BY_KET) == TransformationType::SORT_BY_KET &&
                   ket_of_states[i] != last_ket) {
            blocks.start.push_back(i);
        }

        last_quantum_number_f = quantum_number_f_of_states[i];
        last_quantum_number_m = quantum_number_m_of_states[i];
        last_parity = parity_of_states[i];
        last_ket = ket_of_states[i];
    }

    if ((label & TransformationType::SORT_BY_QUANTUM_NUMBER_F) ==
        TransformationType::SORT_BY_QUANTUM_NUMBER_F) {
        blocks.transformation_type.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
    }
    if ((label & TransformationType::SORT_BY_QUANTUM_NUMBER_M) ==
        TransformationType::SORT_BY_QUANTUM_NUMBER_M) {
        blocks.transformation_type.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
    if ((label & TransformationType::SORT_BY_PARITY) == TransformationType::SORT_BY_PARITY) {
        blocks.transformation_type.push_back(TransformationType::SORT_BY_PARITY);
    }
    if ((label & TransformationType::SORT_BY_KET) == TransformationType::SORT_BY_KET) {
        blocks.transformation_type.push_back(TransformationType::SORT_BY_KET);
    }

    return blocks;
}

template <typename Derived>
std::shared_ptr<Derived> Basis<Derived>::transform(const Sorting &transformation) const {
    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    // Apply the transformation
    transformed->coefficients.matrix = transformed->coefficients.matrix * transformation.matrix;
    for (auto t : transformation.transformation_type) {
        transformed->coefficients.transformation_type.push_back(t);
    }

    {
        auto tmp(transformed->quantum_number_f_of_states);
        for (int i = 0; i < transformation.matrix.size(); ++i) {
            transformed->quantum_number_f_of_states[i] = tmp[transformation.matrix.indices()[i]];
        }
    }
    {
        auto tmp(quantum_number_m_of_states);
        for (int i = 0; i < transformation.matrix.size(); ++i) {
            transformed->quantum_number_m_of_states[i] = tmp[transformation.matrix.indices()[i]];
        }
    }
    {
        auto tmp(parity_of_states);
        for (int i = 0; i < transformation.matrix.size(); ++i) {
            transformed->parity_of_states[i] = tmp[transformation.matrix.indices()[i]];
        }
    }
    {
        auto tmp(ket_of_states);
        for (int i = 0; i < transformation.matrix.size(); ++i) {
            transformed->ket_of_states[i] = tmp[transformation.matrix.indices()[i]];
        }
    }

    return transformed;
}

template <typename Derived>
std::shared_ptr<Derived>
Basis<Derived>::transform(const Transformation<scalar_t> &transformation) const {
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
        bool object_transformation_is_sorting = true;
        for (auto t : coefficients.transformation_type) {
            if (t != TransformationType::SORT_BY_KET &&
                t != TransformationType::SORT_BY_QUANTUM_NUMBER_F &&
                t != TransformationType::SORT_BY_QUANTUM_NUMBER_M &&
                t != TransformationType::SORT_BY_PARITY &&
                t != TransformationType::SORT_BY_ENERGY) {
                object_transformation_is_sorting = false;
                break;
            }
        }
        if (!object_transformation_is_sorting) {
            throw std::runtime_error("If the object was transformed by a different transformation "
                                     "than sorting, it can not be rotated.");
        }
    }

    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    // Apply the transformation
    transformed->coefficients.matrix = coefficients.matrix * transformation.matrix;
    for (auto t : transformation.transformation_type) {
        transformed->coefficients.transformation_type.push_back(t);
    }

    Eigen::SparseMatrix<float> probs =
        (transformation.matrix.cwiseAbs2().transpose()).template cast<float>();

    float tolerance = 1e-16;

    {
        auto map =
            Eigen::Map<const Eigen::VectorXf>(transformed->quantum_number_f_of_states.data(),
                                              transformed->quantum_number_f_of_states.size());
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (size_t i = 0; i < transformed->quantum_number_f_of_states.size(); ++i) {
            if (diff[i] < tolerance) {
                transformed->quantum_number_f_of_states[i] = val[i];
            } else {
                transformed->quantum_number_f_of_states[i] = std::numeric_limits<float>::max();
            }
        }
    }

    {
        auto map =
            Eigen::Map<const Eigen::VectorXf>(transformed->quantum_number_m_of_states.data(),
                                              transformed->quantum_number_m_of_states.size());
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (size_t i = 0; i < transformed->quantum_number_m_of_states.size(); ++i) {
            if (diff[i] < tolerance) {
                transformed->quantum_number_m_of_states[i] = val[i];
            } else {
                transformed->quantum_number_m_of_states[i] = std::numeric_limits<float>::max();
            }
        }
    }

    {
        auto map = Eigen::Map<const Eigen::VectorXi>(transformed->parity_of_states.data(),
                                                     transformed->parity_of_states.size())
                       .template cast<float>();
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (size_t i = 0; i < transformed->parity_of_states.size(); ++i) {
            if (diff[i] < tolerance) {
                transformed->parity_of_states[i] = static_cast<int>(val[i]);
            } else {
                transformed->parity_of_states[i] = std::numeric_limits<int>::max();
            }
        }
    }

    {
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

        std::set<int> ket_of_states_set(transformed->ket_of_states.begin(),
                                        transformed->ket_of_states.end());
        if (ket_of_states_set.size() != transformed->ket_of_states.size()) {
            throw std::runtime_error(
                "Failed to establish a unique mapping between the states and the kets.");
        }
    }

    return transformed;
}

// Explicit instantiations
#include "basis/BasisAtom.hpp"
#include "basis/BasisClassicalLight.hpp"

template class Basis<BasisAtom<float>>;
template class Basis<BasisAtom<double>>;
template class Basis<BasisAtom<std::complex<float>>>;
template class Basis<BasisAtom<std::complex<double>>>;
template class Basis<BasisClassicalLight<float>>;
template class Basis<BasisClassicalLight<double>>;
template class Basis<BasisClassicalLight<std::complex<float>>>;
template class Basis<BasisClassicalLight<std::complex<double>>>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include "utils/hash.hpp"
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

#ifndef DOCTEST_CONFIG_DISABLE

// A derived ket class
class KetDerivedCreator;

class KetDerived : public Ket<float> {
    friend class KetDerivedCreator;
    struct Private {};

public:
    KetDerived(Private, float f, float m, int p, int new_property)
        : Ket<float>(0, f, m, p), new_property(new_property) {}
    std::string get_label() const override { return "my_label"; }
    size_t get_id() const override {
        size_t seed = 0;
        hash::hash_combine(seed, this->quantum_number_f);
        hash::hash_combine(seed, this->quantum_number_m);
        hash::hash_combine(seed, this->parity);
        hash::hash_combine(seed, new_property);
        return seed;
    }
    size_t get_id_for_different_quantum_number_m(float new_quantum_number_m) const override {
        size_t seed = 0;
        hash::hash_combine(seed, this->quantum_number_f);
        hash::hash_combine(seed, new_quantum_number_m);
        hash::hash_combine(seed, this->parity);
        hash::hash_combine(seed, new_property);
        return seed;
    }
    int get_new_property() const { return new_property; }

private:
    int new_property;
};

// Classes for creating an instance of the derived ket class
class KetDerivedCreator {
public:
    KetDerivedCreator(float f, float m, int p, int new_property)
        : f(f), m(m), p(p), new_property(new_property) {}
    std::shared_ptr<const KetDerived> create() const {
        return std::make_shared<const KetDerived>(KetDerived::Private(), f, m, p, new_property);
    }

private:
    float f;
    float m;
    int p;
    int new_property;
};

// A derived basis class
class BasisDerived;

template <>
struct traits::BasisTraits<BasisDerived> {
    using scalar_t = float;
    using real_t = float;
    using ket_t = KetDerived;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

class BasisDerived : public Basis<BasisDerived> {
    friend class BasisDerivedCreator;
    struct Private {};

public:
    BasisDerived(Private, ketvec_t &&kets) : Basis<BasisDerived>(std::move(kets)) {}
    using Type = BasisDerived;
    using ketvec_t = typename traits::BasisTraits<BasisDerived>::ketvec_t;
};

// Classes for creating an instance of the derived basis class
class BasisDerivedCreator {
public:
    BasisDerivedCreator() = default;
    std::shared_ptr<const BasisDerived> create() const {
        std::vector<std::shared_ptr<const KetDerived>> kets;
        kets.reserve(3);
        kets.push_back(KetDerivedCreator(0.5, 0.5, 1, 42).create());
        kets.push_back(KetDerivedCreator(0.5, 0.5, -1, 42).create());
        kets.push_back(KetDerivedCreator(0.5, -0.5, -1, 42).create());
        return std::make_shared<const BasisDerived>(BasisDerived::Private(), std::move(kets));
    }
};

DOCTEST_TEST_CASE("constructing a class derived from basis") {

    // Sort the basis by parity and the m quantum number
    auto tmp = BasisDerivedCreator().create();
    auto basis = tmp->transform(tmp->get_sorter(TransformationType::SORT_BY_PARITY |
                                                TransformationType::SORT_BY_QUANTUM_NUMBER_M));
    int parity = std::numeric_limits<int>::lowest();
    float quantum_number_m = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < basis->get_number_of_states(); ++i) {
        DOCTEST_CHECK(basis->get_parity(i) >= parity);
        DOCTEST_CHECK(basis->get_quantum_number_m(i) >= quantum_number_m);
        parity = basis->get_parity(i);
        quantum_number_m = basis->get_quantum_number_m(i);
    }

    // Check that the blocks are correctly determined
    auto blocks = basis->get_blocks(TransformationType::SORT_BY_PARITY |
                                    TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    DOCTEST_CHECK(blocks.start[0] == 0);
    DOCTEST_CHECK(blocks.start[1] == 1);
    DOCTEST_CHECK(blocks.start[2] == 2);

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_new_property() == 42);
    }
}

#endif // DOCTEST_CONFIG_DISABLE
