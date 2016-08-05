#include "Basisnames.h"

BasisnamesOne::BasisnamesOne() {
}
BasisnamesOne BasisnamesOne::fromStates(std::vector<StateOne> names) {
   BasisnamesOne basisnames;
   basisnames.names_ = names;
   basisnames.dim_ = names.size();
   return basisnames;
}
BasisnamesOne BasisnamesOne::fromFirst(const Configuration &config) {
   StateOne startstate;
   config["n1"] >> startstate.n;
   config["l1"] >> startstate.l;
   config["j1"] >> startstate.j;
   config["m1"] >> startstate.m;

   BasisnamesOne basisnames;
   basisnames._constructedFromFirst = true;
   basisnames.configure(config);
   basisnames.build(startstate, config["species1"].str());
   return basisnames;
}
BasisnamesOne BasisnamesOne::fromFirst(std::shared_ptr<const BasisnamesTwo> basis_two) {
   Configuration config = basis_two->getConf();
   StateOne startstate;
   config["n1"] >> startstate.n;
   config["l1"] >> startstate.l;
   config["j1"] >> startstate.j;
   config["m1"] >> startstate.m;

   BasisnamesOne basisnames;
   basisnames._constructedFromFirst = true;
   basisnames.configure(config);
   basisnames.build(startstate, config["species1"].str(), basis_two, 0);
   return basisnames;
}
BasisnamesOne BasisnamesOne::fromSecond(const Configuration &config) {
    StateOne startstate;
    config["n2"] >> startstate.n;
    config["l2"] >> startstate.l;
    config["j2"] >> startstate.j;
    config["m2"] >> startstate.m;

    BasisnamesOne basisnames;
    basisnames._constructedFromFirst = false;
    basisnames.configure(config);
    basisnames.build(startstate, config["species2"].str());
    return basisnames;
}
BasisnamesOne BasisnamesOne::fromSecond(std::shared_ptr<const BasisnamesTwo> basis_two) {
    Configuration config = basis_two->getConf();
    StateOne startstate;
    config["n2"] >> startstate.n;
    config["l2"] >> startstate.l;
    config["j2"] >> startstate.j;
    config["m2"] >> startstate.m;

    BasisnamesOne basisnames;
    basisnames._constructedFromFirst = false;
    basisnames.configure(config);
    basisnames.build(startstate, config["species2"].str(),basis_two, 1);
    return basisnames;
}
BasisnamesOne BasisnamesOne::fromBoth(const Configuration &config) {
    StateTwo startstate;
    config["n1"] >> startstate.n[0];
    config["l1"] >> startstate.l[0];
    config["j1"] >> startstate.j[0];
    config["m1"] >> startstate.m[0];
    config["n2"] >> startstate.n[1];
    config["l2"] >> startstate.l[1];
    config["j2"] >> startstate.j[1];
    config["m2"] >> startstate.m[1];

    if (config["species1"].str() != config["species2"].str()) {
        std::cout << "BasisnamesOne::fromBoth can only be used if both atoms are of the same species." << std::endl;
        abort();
    }

    BasisnamesOne basisnames;
    basisnames._constructedFromFirst = false;
    basisnames.configure(config);
    basisnames.build(startstate.order(), config["species1"].str());
    return basisnames;

    /*if ((startstate.n[0] == startstate.n[1]) &&
        (startstate.l[0] == startstate.l[1]) &&
        (startstate.j[0] == startstate.j[1]) &&
        (startstate.m[0] == startstate.m[1])) {
        BasisnamesOne basisnames;
        basisnames._constructedFromFirst = false;
        basisnames.configure(config);
        basisnames.build(startstate.first(), config["species1"].str());
        return basisnames;
    } else {

    }*/ // TODO
}
void BasisnamesOne::build(StateTwo startstate, std::string species) {
    states_initial.push_back(startstate.first()); // TODO correct for idx
    states_initial.push_back(startstate.second()); // TODO correct for idx

    conf["species1"] << species;
    conf["n1"] << startstate.n[0];
    conf["l1"] << startstate.l[0];
    conf["j1"] << startstate.j[0];
    conf["m1"] << startstate.m[0];
    conf["n2"] << startstate.n[1];
    conf["l2"] << startstate.l[1];
    conf["j2"] << startstate.j[1];
    conf["m2"] << startstate.m[1];


    std::unordered_set<StateOne> names_set;

    idx_t idx = 0;

    if (delta_l < 0) delta_l = std::fmax(startstate.l[0],startstate.l[1]) + std::fmax(startstate.n[0],startstate.n[1]) + delta_n - 1;
    if (delta_j < 0) delta_j = std::fmax(startstate.j[0],startstate.j[1]) + std::fmax(startstate.n[0],startstate.n[1]) + delta_n - 0.5;
    if (delta_m < 0) delta_m = std::fmax(startstate.m[0],startstate.m[1]) + std::fmax(startstate.n[0],startstate.n[1]) + delta_n - 0.5;

    // loop over quantum numbers of startstate1
    for (int n = std::fmax(0, startstate.n[0] - delta_n); n <= startstate.n[0] + delta_n; ++n) {
        for (int l = std::fmax(0, startstate.l[0] - delta_l); l <= fmin(n-1,startstate.l[0] + delta_l); ++l) {
            for (float j = std::fmax(fabs(l - startstate.s[0]), startstate.j[0] - delta_j); j <= fmin(l + startstate.s[0], startstate.j[0] + delta_j); ++j) {
                for (float m = std::fmax(-j, startstate.m[0] - delta_m); m <= fmin(j, startstate.m[0] + delta_m); ++m) {
                    auto result = names_set.insert(StateOne(idx,n,l,startstate.s[0],j,m));
                    if (result.second) idx++;
                }
            }
        }
    }

    // loop over quantum numbers of startstate2
    for (int n = std::fmax(0, startstate.n[1] - delta_n); n <= startstate.n[1] + delta_n; ++n) {
        for (int l = std::fmax(0, startstate.l[1] - delta_l); l <= fmin(n-1,startstate.l[1] + delta_l); ++l) {
            for (float j = std::fmax(fabs(l - startstate.s[1]), startstate.j[1] - delta_j); j <= fmin(l + startstate.s[1], startstate.j[1] + delta_j); ++j) {
                for (float m = std::fmax(-j, startstate.m[1] - delta_m); m <= fmin(j, startstate.m[1] + delta_m); ++m) {
                    auto result = names_set.insert(StateOne(idx,n,l,startstate.s[1],j,m));
                    if (result.second) idx++;
                }
            }
        }
    }

    std::set<StateOne> names_ordered(names_set.begin(), names_set.end());
    names_ = std::vector<StateOne>(names_ordered.begin(), names_ordered.end());

    dim_ = idx;
}
void BasisnamesOne::build(StateOne startstate, std::string species) {
    states_initial.push_back(startstate); // TODO correct for idx

    conf["species1"] << species;
    conf["n1"] << startstate.n;
    conf["l1"] << startstate.l;
    conf["j1"] << startstate.j;
    conf["m1"] << startstate.m;
    conf["n2"] << "";
    conf["l2"] << "";
    conf["j2"] << "";
    conf["m2"] << "";

    idx_t idx = 0;

    if (delta_l < 0) delta_l = startstate.l + startstate.n + delta_n - 1;
    if (delta_j < 0) delta_j = startstate.j + startstate.n + delta_n - 0.5;
    if (delta_m < 0) delta_m = startstate.m + startstate.n + delta_n - 0.5;

    // loop over quantum numbers
    for (int n = std::fmax(0, startstate.n - delta_n); n <= startstate.n + delta_n; ++n) {
        for (int l = std::fmax(0, startstate.l - delta_l); l <= fmin(n-1,startstate.l + delta_l); ++l) {
            for (float j = std::fmax(fabs(l - startstate.s), startstate.j - delta_j); j <= fmin(l + startstate.s, startstate.j + delta_j); ++j) {
                for (float m = std::fmax(-j, startstate.m - delta_m); m <= fmin(j, startstate.m + delta_m); ++m) { // TODO
                    names_.push_back(StateOne(idx++,n,l,startstate.s,j,m));
                }
            }
        }
    }

    dim_ = idx;
}
void BasisnamesOne::build(StateOne startstate, std::string species, std::shared_ptr<const BasisnamesTwo> basis_two, int i) {
    states_initial.push_back(startstate); // TODO correct for idx

    conf["species1"] << species;
    conf["n1"] << startstate.n;
    conf["l1"] << startstate.l;
    conf["j1"] << startstate.j;
    conf["m1"] << startstate.m;
    conf["n2"] << "";
    conf["l2"] << "";
    conf["j2"] << "";
    conf["m2"] << "";

    std::unordered_set<StateOne> names_set;

    idx_t idx = 0;

    // loop over quantum numbers
    for (auto state : *basis_two) {
        auto result = names_set.insert(StateOne(idx,state.n[i],state.l[i],0.5,state.j[i],state.m[i]));
        if (result.second) idx++;
    }

    std::set<StateOne> names_ordered(names_set.begin(), names_set.end());
    names_ = std::vector<StateOne>(names_ordered.begin(), names_ordered.end());

    dim_ = idx;
}
const std::vector<StateOne>& BasisnamesOne::initial() const {
    return states_initial;
}
void BasisnamesOne::removeUnnecessaryStates(const std::vector<bool> &is_necessary) {
    auto tmp = names_;
    names_.clear();
    names_.reserve(tmp.size());

    // loop over all one-atom states
    idx_t idx = 0;
    for (auto state : tmp) {
        if (is_necessary[state.idx]) {
            state.idx = idx;
            names_.push_back(state);
            // TODO update indices of states_initial
            ++idx;
        }
    }

    dim_ = idx;
    names_.shrink_to_fit();
}


bool BasisnamesOne::constructedFromFirst() {
    return _constructedFromFirst;
}

void BasisnamesOne::save(std::string path) {
    std::ofstream csvfile;
    csvfile.open(path);
    for (const auto &state: *this) {
        csvfile << state.idx << "\t" << state.n << "\t" << state.l << "\t" << state.j << "\t" << state.m <<  std::endl;
    }
    csvfile.close();
}

BasisnamesTwo::BasisnamesTwo(std::shared_ptr<const BasisnamesOne> basis_one1) {
    const Configuration conf1 = basis_one1->getConf();

    if (conf1["n2"].str() == "") {
        std::cout << "BasisnamesTwo can be only constructed from two BasisnamesOne::fromFirst / BasisnamesOne::fromSecond." << std::endl;
        abort();
    }

    configure(conf1);
    conf["combined"] << 1;

    StateTwo startstate;
    conf1["n1"] >> startstate.n[0];
    conf1["l1"] >> startstate.l[0];
    conf1["j1"] >> startstate.j[0];
    conf1["m1"] >> startstate.m[0];
    conf1["n2"] >> startstate.n[1];
    conf1["l2"] >> startstate.l[1];
    conf1["j2"] >> startstate.j[1];
    conf1["m2"] >> startstate.m[1];

    std::array<std::string,2> species({{conf1["species1"].str(),conf1["species1"].str()}}); // TODO : species in state class mit aufnehmen
    build(startstate, species, basis_one1, basis_one1);
}


BasisnamesTwo::BasisnamesTwo(std::shared_ptr<const BasisnamesOne> basis_one1, std::shared_ptr<const BasisnamesOne> basis_one2) {
    const Configuration conf1 = basis_one1->getConf();
    const Configuration conf2 = basis_one2->getConf();

    if (conf1["n2"].str() != "" || conf2["n2"].str() != "") {
        std::cout << "BasisnamesTwo can be only constructed from one single BasisnamesOne::fromBoth." << std::endl;
        abort();
    }

    configure(conf1);
    conf["combined"] << 0;

    StateTwo startstate;
    conf1["n1"] >> startstate.n[0];
    conf1["l1"] >> startstate.l[0];
    conf1["j1"] >> startstate.j[0];
    conf1["m1"] >> startstate.m[0];
    conf2["n1"] >> startstate.n[1];
    conf2["l1"] >> startstate.l[1];
    conf2["j1"] >> startstate.j[1];
    conf2["m1"] >> startstate.m[1];


    std::array<std::string,2> species({{conf1["species1"].str(),conf2["species1"].str()}}); // TODO : species in state class mit aufnehmen
    build(startstate, species, basis_one1, basis_one2);
}

const StateTwo& BasisnamesTwo::initial() const {
    return state_initial;
}

void BasisnamesTwo::removeUnnecessaryStates(const std::vector<bool> &is_necessary) {
    auto tmp = names_;
    names_.clear();
    names_.reserve(tmp.size());

    // loop over all two-atom states
    bool state_initial_found = false;
    idx_t idx = 0;
    for (auto state : tmp) {
        if (is_necessary[state.idx]) {
            state.idx = idx;
            names_.push_back(state);

            if (!state_initial_found && state == state_initial) {
                state_initial.idx = idx;
                state_initial_found = true;
            }

            ++idx;
        }
    }

    dim_ = idx;
    names_.shrink_to_fit();
}

void BasisnamesTwo::build(StateTwo startstate, std::array<std::string,2> species, std::shared_ptr<const BasisnamesOne> basis_one1, std::shared_ptr<const BasisnamesOne> basis_one2) {
    state_initial = startstate;

    conf["species1"] << species[0];
    conf["n1"] << startstate.n[0];
    conf["l1"] << startstate.l[0];
    conf["j1"] << startstate.j[0];
    conf["m1"] << startstate.m[0];
    conf["species2"] << species[1];
    conf["n2"] << startstate.n[1];
    conf["l2"] << startstate.l[1];
    conf["j2"] << startstate.j[1];
    conf["m2"] << startstate.m[1];

    size_t size = basis_one1->size()*basis_one2->size();
    names_.reserve(size);

    idx_t idx = 0;

    // loop over single atom states
    bool state_initial_found = false;
    for (const auto &state_1 : *basis_one1) {
        for (const auto &state_2 : *basis_one2) {
            names_.push_back(StateTwo(idx,state_1,state_2));

            if (!state_initial_found && names_.back() == state_initial) {
                state_initial.idx = idx;
                state_initial_found = true;
            }

            ++idx;
        }
    }

    dim_ = idx;
}

void BasisnamesTwo::save(std::string path) {
    std::ofstream csvfile;
    csvfile.open(path);
    for (const auto &state: *this) {
        csvfile << state.idx << "\t" << state.n[0] << "\t" << state.l[0] << "\t" << state.j[0] << "\t" << state.m[0] << "\t" << state.n[1] << "\t" << state.l[1] << "\t" << state.j[1] << "\t" << state.m[1] << std::endl;
    }
    csvfile.close();
}
