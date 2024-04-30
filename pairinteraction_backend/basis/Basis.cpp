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
    : TransformableSortable<typename traits::BasisTraits<Derived>::scalar_t>(TransformBy::IDENTITY,
                                                                             SortBy::KET),
      kets(std::move(kets)) {
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
    coefficients = Eigen::SparseMatrix<scalar_t>(this->kets.size(), this->kets.size());
    coefficients.setIdentity();
}

template <typename Derived>
const Derived &Basis<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
size_t Basis<Derived>::get_number_of_states() const {
    return coefficients.cols();
}

template <typename Derived>
size_t Basis<Derived>::get_number_of_kets() const {
    return coefficients.rows();
}

template <typename Derived>
const typename Basis<Derived>::ketvec_t &Basis<Derived>::get_kets() const {
    return kets;
}

template <typename Derived>
const Eigen::SparseMatrix<typename Basis<Derived>::scalar_t, Eigen::RowMajor> &
Basis<Derived>::get_coefficients() const {
    return coefficients;
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
const typename Basis<Derived>::ket_t &Basis<Derived>::Iterator::operator*() const {
    return *basis.kets[index];
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++index;
    return *this;
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t>
Basis<Derived>::impl_get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    Eigen::SparseMatrix<scalar_t> rotator(coefficients.rows(), coefficients.rows());

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

    rotator.setFromTriplets(entries.begin(), entries.end());
    rotator.makeCompressed();

    return rotator;
}

template <typename Derived>
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
Basis<Derived>::impl_get_sorter(SortBy label) const {
    std::vector<int> perm(coefficients.cols());
    std::iota(perm.begin(), perm.end(), 0);

    if ((label & SortBy::QUANTUM_NUMBER_F) == SortBy::QUANTUM_NUMBER_F) {
        std::stable_sort(perm.begin(), perm.end(), [&](int i, int j) {
            return quantum_number_f_of_states[i] < quantum_number_f_of_states[j];
        });

        if (quantum_number_f_of_states[perm.back()] == std::numeric_limits<float>::max()) {
            throw std::invalid_argument(
                "States cannot be labeled and thus not sorted by the quantum number f.");
        }
    }

    if ((label & SortBy::QUANTUM_NUMBER_M) == SortBy::QUANTUM_NUMBER_M) {
        std::stable_sort(perm.begin(), perm.end(), [&](int i, int j) {
            return quantum_number_m_of_states[i] < quantum_number_m_of_states[j];
        });

        if (quantum_number_m_of_states[perm.back()] == std::numeric_limits<float>::max()) {
            throw std::invalid_argument(
                "States cannot be labeled and thus not sorted by the quantum number m.");
        }
    }

    if ((label & SortBy::PARITY) == SortBy::PARITY) {
        std::stable_sort(perm.begin(), perm.end(),
                         [&](int i, int j) { return parity_of_states[i] < parity_of_states[j]; });

        if (parity_of_states[perm.back()] == std::numeric_limits<int>::max()) {
            throw std::invalid_argument(
                "States cannot be labeled and thus not sorted by the parity.");
        }
    }

    if ((label & SortBy::KET) == SortBy::KET) {
        std::stable_sort(perm.begin(), perm.end(),
                         [&](int i, int j) { return ket_of_states[i] < ket_of_states[j]; });
    }

    if ((label & SortBy::DIAGONAL) == SortBy::DIAGONAL) {
        throw std::invalid_argument("States do not store the energy and thus can not be sorted by "
                                    "the energy. Use an energy operator instead.");
    }

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> sorter;
    sorter.indices() = Eigen::Map<const Eigen::VectorXi>(perm.data(), perm.size());

    return sorter;
}

template <typename Derived>
std::vector<int> Basis<Derived>::impl_get_blocks(SortBy label) const {
    std::vector<int> blocks;
    float last_quantum_number_f = quantum_number_f_of_states[0];
    float last_quantum_number_m = quantum_number_m_of_states[0];
    int last_parity = parity_of_states[0];
    int last_ket = ket_of_states[0];
    blocks.push_back(0);

    for (int i = 0; i < coefficients.cols(); ++i) {
        if ((label & SortBy::QUANTUM_NUMBER_F) == SortBy::QUANTUM_NUMBER_F &&
            quantum_number_f_of_states[i] != last_quantum_number_f) {
            blocks.push_back(i);
        } else if ((label & SortBy::QUANTUM_NUMBER_M) == SortBy::QUANTUM_NUMBER_M &&
                   quantum_number_m_of_states[i] != last_quantum_number_m) {
            blocks.push_back(i);
        } else if ((label & SortBy::PARITY) == SortBy::PARITY &&
                   parity_of_states[i] != last_parity) {
            blocks.push_back(i);
        } else if ((label & SortBy::KET) == SortBy::KET && ket_of_states[i] != last_ket) {
            blocks.push_back(i);
        }

        last_quantum_number_f = quantum_number_f_of_states[i];
        last_quantum_number_m = quantum_number_m_of_states[i];
        last_parity = parity_of_states[i];
        last_ket = ket_of_states[i];
    }

    return blocks;
}

template <typename Derived>
void Basis<Derived>::impl_transform(const Eigen::SparseMatrix<scalar_t> &transformator) {
    coefficients = coefficients * transformator;

    Eigen::SparseMatrix<float> probs =
        (transformator.cwiseAbs2().transpose()).template cast<float>();

    float tolerance = 1e-16;

    {
        auto map = Eigen::Map<const Eigen::VectorXf>(quantum_number_f_of_states.data(),
                                                     quantum_number_f_of_states.size());
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (size_t i = 0; i < quantum_number_f_of_states.size(); ++i) {
            if (diff[i] < tolerance) {
                quantum_number_f_of_states[i] = val[i];
            } else {
                quantum_number_f_of_states[i] = std::numeric_limits<float>::max();
            }
        }
    }

    {
        auto map = Eigen::Map<const Eigen::VectorXf>(quantum_number_m_of_states.data(),
                                                     quantum_number_m_of_states.size());
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (size_t i = 0; i < quantum_number_m_of_states.size(); ++i) {
            if (diff[i] < tolerance) {
                quantum_number_m_of_states[i] = val[i];
            } else {
                quantum_number_m_of_states[i] = std::numeric_limits<float>::max();
            }
        }
    }

    {
        auto map =
            Eigen::Map<const Eigen::VectorXi>(parity_of_states.data(), parity_of_states.size())
                .template cast<float>();
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (size_t i = 0; i < parity_of_states.size(); ++i) {
            if (diff[i] < tolerance) {
                parity_of_states[i] = static_cast<int>(val[i]);
            } else {
                parity_of_states[i] = std::numeric_limits<int>::max();
            }
        }
    }

    {
        std::vector<real_t> map_idx_to_max(coefficients.cols(), 0);
        for (int row = 0; row < coefficients.outerSize(); ++row) {
            for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(
                     coefficients, row);
                 it; ++it) {
                if (std::abs(it.value()) > map_idx_to_max[it.col()]) {
                    map_idx_to_max[it.col()] = std::abs(it.value());
                    ket_of_states[it.col()] = row;
                }
            }
        }

        std::set<int> ket_of_states_set(ket_of_states.begin(), ket_of_states.end());
        if (ket_of_states_set.size() != ket_of_states.size()) {
            throw std::runtime_error(
                "Failed to establish a unique mapping between the states and the kets.");
        }
    }
}

template <typename Derived>
void Basis<Derived>::impl_sort(
    const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) {
    {
        auto tmp(quantum_number_f_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            quantum_number_f_of_states[i] = tmp[sorter.indices()[i]];
        }
    }
    {
        auto tmp(quantum_number_m_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            quantum_number_m_of_states[i] = tmp[sorter.indices()[i]];
        }
    }
    {
        auto tmp(parity_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            parity_of_states[i] = tmp[sorter.indices()[i]];
        }
    }
    {
        auto tmp(ket_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            ket_of_states[i] = tmp[sorter.indices()[i]];
        }
    }
    coefficients = coefficients * sorter;
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
public:
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
    friend class KetDerivedCreator;
    KetDerived(float f, float m, int p, int new_property)
        : Ket<float>(0, f, m, p), new_property(new_property) {}
    int new_property;
};

// Classes for creating an instance of the derived ket class
class KetDerivedCreator {
public:
    KetDerivedCreator(float f, float m, int p, int new_property)
        : f(f), m(m), p(p), new_property(new_property) {}
    KetDerived create() const { return KetDerived(f, m, p, new_property); }

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
public:
    using Type = BasisDerived;
    using ketvec_t = typename traits::BasisTraits<BasisDerived>::ketvec_t;

private:
    friend class BasisDerivedCreator;
    BasisDerived(ketvec_t &&kets) : Basis<BasisDerived>(std::move(kets)) {}
};

// Classes for creating an instance of the derived basis class
class BasisDerivedCreator {
public:
    BasisDerivedCreator() = default;
    BasisDerived create() const {
        std::vector<std::shared_ptr<const KetDerived>> kets;
        kets.reserve(3);
        kets.push_back(
            std::make_shared<const KetDerived>(KetDerivedCreator(0.5, 0.5, 1, 42).create()));
        kets.push_back(
            std::make_shared<const KetDerived>(KetDerivedCreator(0.5, 0.5, -1, 42).create()));
        kets.push_back(
            std::make_shared<const KetDerived>(KetDerivedCreator(0.5, -0.5, -1, 42).create()));
        return BasisDerived(std::move(kets));
    }
};

DOCTEST_TEST_CASE("constructing a class derived from basis") {
    auto basis = BasisDerivedCreator().create();

    // Sort the basis by parity and the m quantum number
    basis.sort(SortBy::PARITY | SortBy::QUANTUM_NUMBER_M);
    int parity = std::numeric_limits<int>::lowest();
    float quantum_number_m = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < basis.get_number_of_states(); ++i) {
        DOCTEST_CHECK(basis.get_parity(i) >= parity);
        DOCTEST_CHECK(basis.get_quantum_number_m(i) >= quantum_number_m);
        parity = basis.get_parity(i);
        quantum_number_m = basis.get_quantum_number_m(i);
    }

    // Check that the blocks are correctly determined
    auto blocks = basis.get_blocks(SortBy::PARITY | SortBy::QUANTUM_NUMBER_M);
    DOCTEST_CHECK(blocks[0] == 0);
    DOCTEST_CHECK(blocks[1] == 1);
    DOCTEST_CHECK(blocks[2] == 2);

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_new_property() == 42);
    }
}

#endif // DOCTEST_CONFIG_DISABLE
