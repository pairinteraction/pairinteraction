#ifndef BASISNAMES_H
#define BASISNAMES_H

#include "dtypes.h"
#include "Iter.h"

#include <vector>
#include <unordered_set>

template<class T> class Basisnames {
public:
    Basisnames() {}

    size_t size() const {
        return names_.size();
    }
    size_t dim() const {
        return dim_;
    }
    T& get (size_t idx) {
        return names_[idx];
    }
    const T& get (size_t idx) const {
        return names_[idx];
    }
    void set (size_t i, const T &v) {
        names_[i] = v;
    }
    Iter<Basisnames, T> begin() const {
        return Iter<Basisnames, T>( this, 0 );
    }
    Iter<Basisnames, T> end() const {
        return Iter<Basisnames, T>( this, names_.size() );
    }

protected:
    std::vector<T> names_;
    size_t dim_;
};


class BasisnamesOne : public Basisnames<StateOne>{
public:
    BasisnamesOne(const StateOne &startstate);
    BasisnamesOne(const StateOne &startstate1, const StateOne &startstate2);
    const std::vector<StateOne>& initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &is_necessary);
private:
    std::vector<StateOne> states_initial;
};


class BasisnamesTwo : public Basisnames<StateTwo>{
public:
    BasisnamesTwo(const BasisnamesOne &basis_one1, const BasisnamesOne &basis_one2);
    void removeUnnecessaryStates(const std::vector<bool> &isNecessary);
};

#endif
