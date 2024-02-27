#ifndef BASIS_HPP
#define BASIS_HPP

#include <Eigen/SparseCore>
#include <complex>
#include <string>
#include <vector>

#include "state/State.hpp"

template <typename T, bool is_complex>
class Basis {
public:
    Basis();
    T get_energy(size_t index) const;
    float get_quantum_number_f(size_t index) const;
    float get_quantum_number_m(size_t index) const;
    int get_parity(size_t index) const;
    std::string get_label(size_t index) const;

    virtual State<T> get_state(size_t index) const = 0;

    class Iterator {
    public:
        Iterator(const Basis<T, is_complex> &basis, size_t index);
        bool operator!=(const Iterator &other) const;
        State<T> operator*() const;
        Iterator &operator++();

    private:
        const Basis<T, is_complex> &basis;
        size_t index;
    };

    Iterator begin() const;
    Iterator end() const;

protected:
    void reserve_coefficients(int n);
    void add_coefficients(T energy, float f, float m, int p, std::string label);
    void assemble_coefficients();
    bool is_assembled;
    bool is_standard_basis;

private:
    using scalar_t = typename std::conditional<is_complex, std::complex<T>, T>::type;
    std::vector<T> energies;
    std::vector<float> quantum_numbers_f;
    std::vector<float> quantum_numbers_m;
    std::vector<int> parities;
    std::vector<std::string> labels;
    Eigen::SparseMatrix<scalar_t> coefficients;
    size_t num_coefficients;
};

#endif // BASIS_HPP
