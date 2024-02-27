#ifndef STATE_HPP
#define STATE_HPP

#include <string>

template <typename T>
class State {
public:
    State(T energy, float f, float m, int p, std::string label);
    T get_energy() const;
    float get_quantum_number_f() const;
    float get_quantum_number_m() const;
    int get_parity() const;
    std::string get_label() const;

private:
    T energy;
    float quantum_number_f;
    float quantum_number_m;
    int parity;
    std::string label;
};

#endif // STATE_HPP
