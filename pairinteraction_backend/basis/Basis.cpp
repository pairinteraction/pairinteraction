#include "basis/Basis.hpp"
#include "utils/euler.hpp"

#include <numeric>
#include <set>

template <typename Derived>
Basis<Derived>::Basis(ketvec_t &&kets)
    : kets(std::move(kets)), is_standard_basis(true), sortation(Label::KET) {
    quantum_number_f_of_states.reserve(this->kets.size());
    quantum_number_m_of_states.reserve(this->kets.size());
    parity_of_states.reserve(this->kets.size());
    for (const auto &ket : this->kets) {
        quantum_number_f_of_states.push_back(ket->get_quantum_number_f());
        quantum_number_m_of_states.push_back(ket->get_quantum_number_m());
        parity_of_states.push_back(ket->get_parity());
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
const typename Basis<Derived>::ket_t &Basis<Derived>::get_ket(size_t index_ket) const {
    return *kets[index_ket];
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
    return basis.get_ket(index);
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++index;
    return *this;
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t>
Basis<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t>
Basis<Derived>::get_rotator(std::array<real_t, 3> to_z_axis,
                            std::array<real_t, 3> to_y_axis) const {
    auto euler_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->get_rotator(euler_angles[0], euler_angles[1], euler_angles[2]);
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter(Label label) const {
    std::vector<int> sorter(coefficients.cols());
    std::iota(sorter.begin(), sorter.end(), 0);

    if ((label & Label::QUANTUM_NUMBER_F) == Label::QUANTUM_NUMBER_F) {
        std::stable_sort(sorter.begin(), sorter.end(), [&](int i, int j) {
            return quantum_number_f_of_states[i] < quantum_number_f_of_states[j];
        });
    }

    if ((label & Label::QUANTUM_NUMBER_M) == Label::QUANTUM_NUMBER_M) {
        std::stable_sort(sorter.begin(), sorter.end(), [&](int i, int j) {
            return quantum_number_m_of_states[i] < quantum_number_m_of_states[j];
        });
    }

    if ((label & Label::PARITY) == Label::PARITY) {
        std::stable_sort(sorter.begin(), sorter.end(),
                         [&](int i, int j) { return parity_of_states[i] < parity_of_states[j]; });
    }

    if ((label & Label::KET) == Label::KET) {
        std::stable_sort(sorter.begin(), sorter.end(),
                         [&](int i, int j) { return ket_of_states[i] < ket_of_states[j]; });
    }

    return sorter;
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_blocks(Label label) const {
    if (label != sortation) {
        throw std::invalid_argument("The basis is not sorted by the requested label.");
    }

    std::vector<int> blocks;
    float last_quantum_number_f = quantum_number_f_of_states[0];
    float last_quantum_number_m = quantum_number_m_of_states[0];
    int last_parity = parity_of_states[0];
    int last_ket = ket_of_states[0];
    blocks.push_back(0);

    for (size_t i = 0; i < coefficients.cols(); ++i) {
        if ((label & Label::QUANTUM_NUMBER_F) == Label::QUANTUM_NUMBER_F &&
            quantum_number_f_of_states[i] != last_quantum_number_f) {
            blocks.push_back(i);
        } else if ((label & Label::QUANTUM_NUMBER_M) == Label::QUANTUM_NUMBER_M &&
                   quantum_number_m_of_states[i] != last_quantum_number_m) {
            blocks.push_back(i);
        } else if ((label & Label::PARITY) == Label::PARITY && parity_of_states[i] != last_parity) {
            blocks.push_back(i);
        } else if ((label & Label::KET) == Label::KET && ket_of_states[i] != last_ket) {
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
void Basis<Derived>::transform(const Eigen::SparseMatrix<scalar_t> &transformator) {
    is_standard_basis = false;
    coefficients = coefficients * transformator;
    sortation = Label::NONE;

    float tolerance = 1e-16;

    Eigen::SparseMatrix<scalar_t> identity(transformator.cols(), transformator.cols());
    identity.setIdentity();
    if ((coefficients.adjoint() * coefficients - identity).norm() > tolerance) {
        throw std::runtime_error("The transformation is not unitary.");
    }

    Eigen::SparseMatrix<float> probs =
        (transformator.cwiseAbs2().transpose()).template cast<float>();

    {
        auto map = Eigen::Map<const Eigen::VectorXf>(quantum_number_f_of_states.data(),
                                                     quantum_number_f_of_states.size());
        Eigen::VectorXf val = probs * map;
        Eigen::VectorXf sq = probs * map.cwiseAbs2();
        Eigen::VectorXf diff = (val * val - sq).cwiseAbs();

        for (int i = 0; i < quantum_number_f_of_states.size(); ++i) {
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

        for (int i = 0; i < quantum_number_m_of_states.size(); ++i) {
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

        for (int i = 0; i < parity_of_states.size(); ++i) {
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
            for (typename Eigen::SparseMatrix<scalar_t>::InnerIterator it(coefficients, row); it;
                 ++it) {
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
void Basis<Derived>::rotate(real_t alpha, real_t beta, real_t gamma) {
    auto rotator = this->get_rotator(alpha, beta, gamma);
    this->transform(rotator);
}

template <typename Derived>
void Basis<Derived>::rotate(std::array<real_t, 3> to_z_axis, std::array<real_t, 3> to_y_axis) {
    auto rotator = this->get_rotator(to_z_axis, to_y_axis);
    this->transform(rotator);
}

template <typename Derived>
void Basis<Derived>::sort(const std::vector<int> &sorter) {
    {
        auto tmp(quantum_number_f_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            quantum_number_f_of_states[i] = tmp[sorter[i]];
        }
    }
    {
        auto tmp(quantum_number_m_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            quantum_number_m_of_states[i] = tmp[sorter[i]];
        }
    }
    {
        auto tmp(parity_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            parity_of_states[i] = tmp[sorter[i]];
        }
    }
    {
        auto tmp(ket_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            ket_of_states[i] = tmp[sorter[i]];
        }
    }
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
    perm.indices() = Eigen::Map<const Eigen::VectorXi>(sorter.data(), sorter.size());
    coefficients = coefficients * perm;
    sortation = Label::NONE;
}

template <typename Derived>
void Basis<Derived>::sort(Label label) {
    auto sorter = this->get_sorter(label);
    this->sort(sorter);
    sortation = label;
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

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

#ifndef DOCTEST_CONFIG_DISABLE

// A derived ket class
class KetDerived : public Ket<float> {
public:
    int get_new_property() const { return new_property; }

private:
    friend class KetDerivedCreator;
    KetDerived(float f, float m, int p, std::string label, int new_property)
        : Ket<float>(0, f, m, p, label), new_property(new_property) {}
    int new_property;
};

// Classes for creating an instance of the derived ket class
class KetDerivedCreator {
public:
    KetDerivedCreator(float f, float m, int p, std::string label, int new_property)
        : f(f), m(m), p(p), label(label), new_property(new_property) {}
    KetDerived create() const { return KetDerived(f, m, p, label, new_property); }

private:
    float f;
    float m;
    int p;
    std::string label;
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
            std::make_shared<const KetDerived>(KetDerivedCreator(0.5, 0.5, 1, "1s", 42).create()));
        kets.push_back(
            std::make_shared<const KetDerived>(KetDerivedCreator(0.5, 0.5, -1, "2s", 42).create()));
        kets.push_back(std::make_shared<const KetDerived>(
            KetDerivedCreator(0.5, -0.5, -1, "3s", 42).create()));
        return BasisDerived(std::move(kets));
    }
};

DOCTEST_TEST_CASE("constructing a class derived from basis") {
    auto basis = BasisDerivedCreator().create();

    // Sort the basis by parity and the m quantum number
    basis.sort(BasisDerived::Label::PARITY | BasisDerived::Label::QUANTUM_NUMBER_M);
    int parity = std::numeric_limits<int>::lowest();
    float quantum_number_m = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < basis.get_number_of_states(); ++i) {
        DOCTEST_CHECK(basis.get_parity(i) >= parity);
        DOCTEST_CHECK(basis.get_quantum_number_m(i) >= quantum_number_m);
        parity = basis.get_parity(i);
        quantum_number_m = basis.get_quantum_number_m(i);
    }

    // Check that the blocks are correctly determined
    auto blocks =
        basis.get_blocks(BasisDerived::Label::PARITY | BasisDerived::Label::QUANTUM_NUMBER_M);
    DOCTEST_CHECK(blocks[0] == 0);
    DOCTEST_CHECK(blocks[1] == 1);
    DOCTEST_CHECK(blocks[2] == 2);

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_new_property() == 42);
    }
}

#endif // DOCTEST_CONFIG_DISABLE
