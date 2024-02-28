// #pragma once

// #include "Basis.hpp"

// #include <limits>

// template <typename T, bool is_complex>
// class BasisAtom : public Basis<T, is_complex> {
// public:
//     BasisAtom();
//     void restrict_energy(T min, T max);
//     void restrict_quantum_number_f(float min, float max);
//     void restrict_quantum_number_m(float min, float max);
//     void restrict_parity(int parity);
//     void restrict_quantum_number_n(int n);
//     void restrict_quantum_number_nu(T min, T max);
//     void restrict_quantum_number_l(T min, T max);
//     void restrict_quantum_number_s(T min, T max);
//     void restrict_quantum_number_j(T min, T max);

// private:
//     void ensure_assembled_kets() override;
//     T min_energy{std::numeric_limits<T>::lowest()};
//     T max_energy{std::numeric_limits<T>::max()};
//     float min_quantum_number_f{std::numeric_limits<float>::lowest()};
//     float max_quantum_number_f{std::numeric_limits<float>::max()};
//     float min_quantum_number_m{std::numeric_limits<float>::lowest()};
//     float max_quantum_number_m{std::numeric_limits<float>::max()};
//     int parity{std::numeric_limits<int>::max()};
//     int quantum_number_n{std::numeric_limits<int>::max()};
//     T min_quantum_number_nu{std::numeric_limits<T>::lowest()};
//     T max_quantum_number_nu{std::numeric_limits<T>::max()};
//     T min_quantum_number_l{std::numeric_limits<T>::lowest()};
//     T max_quantum_number_l{std::numeric_limits<T>::max()};
//     T min_quantum_number_s{std::numeric_limits<T>::lowest()};
//     T max_quantum_number_s{std::numeric_limits<T>::max()};
//     T min_quantum_number_j{std::numeric_limits<T>::lowest()};
//     T max_quantum_number_j{std::numeric_limits<T>::max()};
// };

// extern template class BasisAtom<float, false>;
// extern template class BasisAtom<double, false>;
// extern template class BasisAtom<float, true>;
// extern template class BasisAtom<double, true>;
