#ifndef BASIS_HPP
#define BASIS_HPP

#include <Eigen/SparseCore>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "ket/Ket.hpp"

template <typename T, bool is_complex>
class Basis {
public:
    Basis();
    const Ket<T> &get_ket(size_t index);
    T get_energy(size_t index);
    float get_quantum_number_f(size_t index);
    float get_quantum_number_m(size_t index);
    int get_parity(size_t index);
    std::string get_label(size_t index);

    class Iterator {
    public:
        Iterator(const Basis<T, is_complex> &basis, size_t index);
        bool operator!=(const Iterator &other) const;
        const Ket<T> &operator*() const;
        Iterator &operator++();

    private:
        const Basis<T, is_complex> &basis;
        size_t index;
    };

    Iterator begin();
    Iterator end();

protected:
    virtual void ensure_assembled_kets() = 0;
    void ensure_assembled();
    void ensure_not_assembled() const;
    void ensure_standard_basis() const;
    std::vector<std::shared_ptr<Ket<T>>> kets; // TODO const Ket?

private:
    using scalar_t = typename std::conditional<is_complex, std::complex<T>, T>::type;
    std::vector<T> energies;
    std::vector<float> quantum_numbers_f;
    std::vector<float> quantum_numbers_m;
    std::vector<int> parities;
    std::vector<std::string> labels;
    Eigen::SparseMatrix<scalar_t> coefficients;
    bool is_assembled;
    bool is_standard_basis;
};

#endif // BASIS_HPP
