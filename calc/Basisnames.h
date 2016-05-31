#ifndef BASISNAMES_H
#define BASISNAMES_H

#include "dtypes.h"
#include "Iter.h"
#include "State.h"
#include "ConfParser.hpp"

#include <vector>
#include <set>
#include <unordered_set>
#include <memory>
#include <string>
#include <fstream>

template<class T> class Basisnames {
public:
    Basisnames() {
    }
    void configure(const Configuration &config) {
        conf["deltaNSingle"] = config["deltaNSingle"];
        conf["deltaLSingle"] = config["deltaLSingle"];
        conf["deltaJSingle"] = config["deltaJSingle"];
        conf["deltaMSingle"] = config["deltaMSingle"];
        conf["deltaNSingle"] >> delta_n;
        conf["deltaLSingle"] >> delta_l;
        conf["deltaJSingle"] >> delta_j;
        conf["deltaMSingle"] >> delta_m;
    }
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
    const Configuration& getConf() const { // TODO in Configurable Klasse auslagern, von der geerbt werrden soll
        return conf;
    }

protected:
    int delta_n, delta_l, delta_j, delta_m;
    Configuration conf;
    std::vector<T> names_;
    size_t dim_;

};

class BasisnamesTwo;


class BasisnamesOne : public Basisnames<StateOne>{
public:
    /*BasisnamesOne(const Configuration &config, const StateOne &startstate);
    BasisnamesOne(const Configuration &config, const StateTwo &startstate);*/
    static BasisnamesOne fromStates(std::vector<StateOne> names); // TODO
    static BasisnamesOne fromFirst(const Configuration &config);
    static BasisnamesOne fromFirst(std::shared_ptr<const BasisnamesTwo> basis_two);
    static BasisnamesOne fromSecond(const Configuration &config);
    static BasisnamesOne fromSecond(std::shared_ptr<const BasisnamesTwo> basis_two);
    static BasisnamesOne fromBoth(const Configuration &config);
    const std::vector<StateOne>& initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &is_necessary);
    bool constructedFromFirst();
    void save(std::string path);
private:
    BasisnamesOne();
    void build(StateOne startstate, std::string species);
    void build(StateOne startstate, std::string species, std::shared_ptr<const BasisnamesTwo> basis_two, int i);
    void build(StateTwo startstate, std::string species);
    std::vector<StateOne> states_initial;
    bool _constructedFromFirst;
};




class BasisnamesTwo : public Basisnames<StateTwo>{
public:
    BasisnamesTwo(std::shared_ptr<const BasisnamesOne> basis_one1, std::shared_ptr<const BasisnamesOne> basis_one2);
    BasisnamesTwo(std::shared_ptr<const BasisnamesOne> basis_one1);
    const StateTwo& initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &isNecessary);
    void save(std::string path);
protected:
    void build(StateTwo startstate, std::array<std::string,2> species, std::shared_ptr<const BasisnamesOne> basis_one1, std::shared_ptr<const BasisnamesOne> basis_one2);
private:
    StateTwo state_initial;
};

#endif
