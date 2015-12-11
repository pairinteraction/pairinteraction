#include "Basisnames.h"

BasisnamesOne::BasisnamesOne(const StateOne &startstate) {
    states_initial.push_back(startstate);

    size_t size = 4; // TODO
    names_.reserve(size);

    idx_t idx = 0;

    int delta_n = 4; // 4 TODO
    int delta_l = 20; // 100 TODO
    int delta_j = 100; // TODO
    int delta_m = 1; // TODO

    // loop over quantum numbers
    for (int n = fmax(0, startstate.n - delta_n); n <= startstate.n + delta_n; ++n) {
        for (int l = fmax(0, startstate.l - delta_l); l <= fmin(n-1,startstate.l + delta_l); ++l) {
            for (float j = fmax(fabs(l - startstate.s), startstate.j - delta_j); j <= fmin(l + startstate.s, startstate.j + delta_j); ++j) {
                //for (float m = fmax(-j, startstate.m - delta_m); m <= fmin(j, startstate.m + delta_m); ++m) {
                for (float m = fmax(-j, -1.5); m <= fmin(j, 1.5); ++m) {
                    names_.push_back(StateOne(idx++,n,l,startstate.s,j,m));
                }
            }
        }
    }

    dim_ = idx;
    std::cout << dim_ << std::endl;
}

BasisnamesOne::BasisnamesOne(const StateOne &startstate1, const StateOne &startstate2) {
    states_initial.push_back(startstate1);
    states_initial.push_back(startstate2);

    size_t size = 4; // TODO
    std::unordered_set<StateOne> names_set; // TODO auf das sortierte set wechseln
    names_set.reserve(size);

    idx_t idx = 0;

    int delta_n = 1; // TODO
    int delta_l = 1; // TODO
    int delta_j = 1; // TODO
    int delta_m = 1; // TODO

    // loop over quantum numbers of startstate1
    for (int n = fmax(0, startstate1.n - delta_n); n <= startstate1.n + delta_n; ++n) {
        for (int l = fmax(0, startstate1.l - delta_l); l <= fmin(n-1,startstate1.l + delta_l); ++l) {
            for (float j = fmax(fabs(l - startstate1.s), startstate1.j - delta_j); j <= fmin(l + startstate1.s, startstate1.j + delta_j); ++j) {
                for (float m = fmax(-j, startstate1.m - delta_m); m <= fmin(j, startstate1.m + delta_m); ++m) {
                    auto result = names_set.insert(StateOne(idx,n,l,startstate1.s,j,m));
                    if (result.second) idx++;
                }
            }
        }
    }

    // loop over quantum numbers of startstate2
    for (int n = fmax(0, startstate2.n - delta_n); n <= startstate2.n + delta_n; ++n) {
        for (int l = fmax(0, startstate2.l - delta_l); l <= fmin(n-1,startstate2.l + delta_l); ++l) {
            for (float j = fmax(fabs(l - startstate2.s), startstate2.j - delta_j); j <= fmin(l + startstate2.s, startstate2.j + delta_j); ++j) {
                for (float m = fmax(-j, startstate2.m - delta_m); m <= fmin(j, startstate2.m + delta_m); ++m) {
                    auto result = names_set.insert(StateOne(idx,n,l,startstate2.s,j,m));
                    if (result.second) idx++;
                }
            }
        }
    }

    dim_ = idx;
    names_ = std::vector<StateOne>(names_set.begin(), names_set.end());
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
