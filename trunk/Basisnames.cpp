#include "Basisnames.h"

BasisnamesOne::BasisnamesOne() {
}
BasisnamesOne BasisnamesOne::fromFirst(const Configuration &config) {
   StateOne startstate;
   config["n1"] >> startstate.n;
   config["l1"] >> startstate.l;
   config["j1"] >> startstate.j;
   config["m1"] >> startstate.m;

   BasisnamesOne basisnames;
   basisnames.configure(config);
   basisnames.build(startstate, config["species1"].str());
   return basisnames;
}
BasisnamesOne BasisnamesOne::fromSecond(const Configuration &config) {
    StateOne startstate;
    config["n2"] >> startstate.n;
    config["l2"] >> startstate.l;
    config["j2"] >> startstate.j;
    config["m2"] >> startstate.m;

    BasisnamesOne basisnames;
    basisnames.configure(config);
    basisnames.build(startstate, config["species2"].str());
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

    if ((startstate.n[0] == startstate.n[1]) &&
        (startstate.l[0] == startstate.l[1]) &&
        (startstate.j[0] == startstate.j[1]) &&
        (startstate.m[0] == startstate.m[1])) {
        BasisnamesOne basisnames;
        basisnames.configure(config);
        basisnames.build(startstate.first(), config["species1"].str());
        return basisnames;
    } else {
        BasisnamesOne basisnames;
        basisnames.configure(config);
        basisnames.build(startstate.order(), config["species1"].str());
        return basisnames;
    }
}
void BasisnamesOne::build(StateTwo startstate, std::string species) {
    states_initial.push_back(startstate.first());
    states_initial.push_back(startstate.second());

    conf["species"] << species;
    conf["n1"] << startstate.n[0];
    conf["l1"] << startstate.l[0];
    conf["j1"] << startstate.j[0];
    conf["m1"] << startstate.m[0];
    conf["n2"] << startstate.n[1];
    conf["l2"] << startstate.l[1];
    conf["j2"] << startstate.j[1];
    conf["m2"] << startstate.m[1];

    std::unordered_set<StateOne> names_set; // TODO auf das sortierte set wechseln

    idx_t idx = 0;

    // loop over quantum numbers of startstate1
    for (int n = fmax(0, startstate.n[0] - delta_n); n <= startstate.n[0] + delta_n; ++n) {
        for (int l = fmax(0, startstate.l[0] - delta_l); l <= fmin(n-1,startstate.l[0] + delta_l); ++l) {
            for (float j = fmax(fabs(l - startstate.s[0]), startstate.j[0] - delta_j); j <= fmin(l + startstate.s[0], startstate.j[0] + delta_j); ++j) {
                for (float m = fmax(-j, startstate.m[0] - delta_m); m <= fmin(j, startstate.m[0] + delta_m); ++m) {
                    auto result = names_set.insert(StateOne(idx,n,l,startstate.s[0],j,m));
                    if (result.second) idx++;
                }
            }
        }
    }

    // loop over quantum numbers of startstate2
    for (int n = fmax(0, startstate.n[1] - delta_n); n <= startstate.n[1] + delta_n; ++n) {
        for (int l = fmax(0, startstate.l[1] - delta_l); l <= fmin(n-1,startstate.l[1] + delta_l); ++l) {
            for (float j = fmax(fabs(l - startstate.s[1]), startstate.j[1] - delta_j); j <= fmin(l + startstate.s[1], startstate.j[1] + delta_j); ++j) {
                for (float m = fmax(-j, startstate.m[1] - delta_m); m <= fmin(j, startstate.m[1] + delta_m); ++m) {
                    auto result = names_set.insert(StateOne(idx,n,l,startstate.s[1],j,m));
                    if (result.second) idx++;
                }
            }
        }
    }

    names_ = std::vector<StateOne>(names_set.begin(), names_set.end());

    dim_ = idx;

    std::cout << dim_ << std::endl;
}
void BasisnamesOne::build(StateOne startstate, std::string species) {
    states_initial.push_back(startstate);

    conf["species"] << species;
    conf["n1"] << startstate.n;
    conf["l1"] << startstate.l;
    conf["j1"] << startstate.j;
    conf["m1"] << startstate.m;
    conf["n2"] << "";
    conf["l2"] << "";
    conf["j2"] << "";
    conf["m2"] << "";

    idx_t idx = 0;

    // loop over quantum numbers
    for (int n = fmax(0, startstate.n - delta_n); n <= startstate.n + delta_n; ++n) {
        for (int l = fmax(0, startstate.l - delta_l); l <= fmin(n-1,startstate.l + delta_l); ++l) {
            for (float j = fmax(fabs(l - startstate.s), startstate.j - delta_j); j <= fmin(l + startstate.s, startstate.j + delta_j); ++j) {
                for (float m = fmax(-j, startstate.m - delta_m); m <= fmin(j, startstate.m + delta_m); ++m) { // TODO
                    names_.push_back(StateOne(idx++,n,l,startstate.s,j,m));
                }
            }
        }
    }

    dim_ = idx;
}
const std::vector<StateOne>& BasisnamesOne::initial() const {
    return states_initial;
}
void BasisnamesOne::removeUnnecessaryStates(const std::vector<bool> &is_necessary) {
    auto tmp = names_;
    names_.clear();
    names_.reserve(tmp.size());

    // loop over all two-atom states
    idx_t idx = 0;
    for (auto state : tmp) {
        if (is_necessary[state.idx]) {
            state.idx = idx++;
            names_.push_back(state);
        }
    }

    dim_ = idx;
    names_.shrink_to_fit();
}





BasisnamesTwo::BasisnamesTwo(const BasisnamesOne &basis_one1, const BasisnamesOne &basis_one2) {
    const Configuration conf1 = basis_one2.getConf();
    const Configuration conf2 = basis_one2.getConf();
    if (conf1 == conf2) {
        conf = conf1;
    } else {
        conf = conf1;
        conf["n1"] = conf1["n"];
        conf["l1"] = conf1["l"];
        conf["j1"] = conf1["j"];
        conf["m1"] = conf1["m"];
        conf["n2"] = conf2["n"];
        conf["l2"] = conf2["l"];
        conf["j2"] = conf2["j"];
        conf["m2"] = conf2["m"];
        // conf.erase("n");
        // conf.erase("l");
        // conf.erase("j");
        // conf.erase("m");
    }

    configure(basis_one1.getConf());

    size_t size = basis_one1.size()*basis_one2.size();
    names_.reserve(size);

    idx_t idx = 0;

    // loop over single atom states
    for (const auto &state_1 : basis_one1) {
        for (const auto &state_2 : basis_one2) {
            names_.push_back(StateTwo(idx++,state_1,state_2));
        }
    }
}

void BasisnamesTwo::removeUnnecessaryStates(const std::vector<bool> &isNecessary) {
    auto tmp = names_;
    names_.clear();
    names_.reserve(tmp.size());

    // loop over all two-atom states
    for (auto state : tmp) {
        if (isNecessary[state.idx]) {
            names_.push_back(state);
        }
    }

    names_.shrink_to_fit();
}
