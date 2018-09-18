/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "PerturbativeInteraction.h"

#include <unordered_set>

PerturbativeInteraction::PerturbativeInteraction(MatrixElementCache &cache)
    : cache(cache), bfield(0) { // TODO remove this constructor
    initializeAngleTerms(0);
}

PerturbativeInteraction::PerturbativeInteraction(double angle, MatrixElementCache &cache)
    : cache(cache), bfield(0) {
    initializeAngleTerms(angle);
}

PerturbativeInteraction::PerturbativeInteraction(double angle, double weak_bfield_along_z,
                                                 MatrixElementCache &cache)
    : cache(cache), bfield(weak_bfield_along_z) {
    initializeAngleTerms(angle);
}

double PerturbativeInteraction::getC6(const StateTwo &state, double deltaN) {
    double C6 = 0;

    for (int n0 = state.getN(0) - deltaN; n0 <= state.getN(0) + deltaN; ++n0) {

        for (int n1 = state.getN(1) - deltaN; n1 <= state.getN(1) + deltaN; ++n1) {

            std::vector<int> array_l0;
            if (state.getL(0) > 0) {
                array_l0.push_back(state.getL(0) - 1);
            }
            if (state.getL(0) < n0 - 1) {
                array_l0.push_back(state.getL(0) + 1);
            }
            for (int l0 : array_l0) {

                std::vector<int> array_l1;
                if (state.getL(1) > 0) {
                    array_l1.push_back(state.getL(1) - 1);
                }
                if (state.getL(1) < n1 - 1) {
                    array_l1.push_back(state.getL(1) + 1);
                }
                for (int l1 : array_l1) {

                    std::set<float> array_j0;
                    if (std::abs(std::abs(l0 - state.getS(0)) - state.getJ(0)) < 2) {
                        array_j0.insert(std::abs(l0 - state.getS(0)));
                    }
                    if (std::abs(l0 + state.getS(0) - state.getJ(0)) < 2) {
                        array_j0.insert(l0 + state.getS(0));
                    }
                    for (float j0 : array_j0) {

                        std::set<float> array_j1;
                        if (std::abs(std::abs(l1 - state.getS(1)) - state.getJ(1)) < 2) {
                            array_j1.insert(std::abs(l1 - state.getS(1)));
                        }
                        if (std::abs(l1 + state.getS(1) - state.getJ(1)) < 2) {
                            array_j1.insert(l1 + state.getS(1));
                        }
                        for (float j1 : array_j1) {

                            for (auto &q : array_q) {
                                // q = final.m - initial.m

                                float m0 = state.getM(0) + q[0];
                                if (std::abs(m0) > j0) {
                                    continue;
                                }

                                float m1 = state.getM(1) + q[1];
                                if (std::abs(m1) > j1) {
                                    continue;
                                }

                                StateTwo state_virtual =
                                    StateTwo(state.getSpecies(), {{n0, n1}}, {{l0, l1}}, {{j0, j1}},
                                             {{m0, m1}});

                                double energydiff = state.getEnergy() - state_virtual.getEnergy();
                                if (bfield != 0) {
                                    energydiff += -bfield *
                                        (cache.getMagneticDipole(state.getFirstState(),
                                                                 state.getFirstState()) +
                                         cache.getMagneticDipole(state.getSecondState(),
                                                                 state.getSecondState()) -
                                         cache.getMagneticDipole(state_virtual.getFirstState(),
                                                                 state_virtual.getFirstState()) -
                                         cache.getMagneticDipole(state_virtual.getSecondState(),
                                                                 state_virtual.getSecondState()));
                                }

                                C6 +=
                                    std::pow(
                                        coulombs_constant *
                                            array_angle_term[3 * (q[0] + 1) + (q[1] + 1)] *
                                            cache.getElectricDipole(state_virtual.getFirstState(),
                                                                    state.getFirstState()) *
                                            cache.getElectricDipole(state_virtual.getSecondState(),
                                                                    state.getSecondState()),
                                        2) /
                                    energydiff; // getMultipole(final, inital, 1)
                            }
                        }
                    }
                }
            }
        }
    }

    return C6;
}

eigen_dense_double_t PerturbativeInteraction::getC6(const std::vector<StateTwo> &states,
                                                    double deltaN) {
    eigen_dense_double_t C6_matrix = eigen_dense_double_t::Zero(states.size(), states.size());

    std::unordered_set<StateTwo> set_states(states.begin(), states.end());

    for (size_t idx_row = 0; idx_row < states.size(); ++idx_row) {
        auto &state_row = states[idx_row];

        for (size_t idx_col = idx_row; idx_col < states.size(); ++idx_col) {
            auto &state_col = states[idx_col];

            if (state_row.getS(0) != state_col.getS(0) || state_row.getS(1) != state_col.getS(1)) {
                continue;
            }
            auto s = state_col.getS();

            if (state_row.getSpecies(0) != state_col.getSpecies(0) ||
                state_row.getSpecies(1) != state_col.getSpecies(1)) {
                continue;
            }
            auto species = state_col.getSpecies();

            int n0_min = std::min(state_row.getN(0), state_col.getN(0));
            int n0_max = std::max(state_row.getN(0), state_col.getN(0));
            int n1_min = std::min(state_row.getN(1), state_col.getN(1));
            int n1_max = std::max(state_row.getN(1), state_col.getN(1));

            double C6 = 0;

            for (int n0 = n0_min - deltaN; n0 <= n0_max + deltaN; ++n0) {

                for (int n1 = n1_min - deltaN; n1 <= n1_max + deltaN; ++n1) {

                    std::vector<int> array_l0;
                    if (state_row.getL(0) == state_col.getL(0)) {
                        if (state_col.getL(0) < n0 - 1) {
                            array_l0.push_back(state_col.getL(0) + 1);
                        }
                        if (state_col.getL(0) > 0) {
                            array_l0.push_back(state_col.getL(0) - 1);
                        }
                    } else if (state_row.getL(0) + 2 == state_col.getL(0) &&
                               state_col.getL(0) < n0 - 1) {
                        array_l0.push_back(state_col.getL(0) + 1);
                    } else if (state_row.getL(0) - 2 == state_col.getL(0) &&
                               state_col.getL(0) > 0) {
                        array_l0.push_back(state_col.getL(0) - 1);
                    }
                    for (int l0 : array_l0) {

                        std::vector<int> array_l1;
                        if (state_row.getL(1) == state_col.getL(1)) {
                            if (state_col.getL(1) < n1 - 1) {
                                array_l1.push_back(state_col.getL(1) + 1);
                            }
                            if (state_col.getL(1) > 0) {
                                array_l1.push_back(state_col.getL(1) - 1);
                            }
                        } else if (state_row.getL(1) + 2 == state_col.getL(1) &&
                                   state_col.getL(1) < n1 - 1) {
                            array_l1.push_back(state_col.getL(1) + 1);
                        } else if (state_row.getL(1) - 2 == state_col.getL(1) &&
                                   state_col.getL(1) > 0) {
                            array_l1.push_back(state_col.getL(1) - 1);
                        }
                        for (int l1 : array_l1) {

                            std::set<float> array_j0;
                            if (std::abs(std::abs(l0 - s[0]) - state_col.getJ(0)) < 2 &&
                                std::abs(std::abs(l0 - s[0]) - state_row.getJ(0)) < 2) {
                                array_j0.insert(std::abs(l0 - s[0]));
                            }
                            if (std::abs(l0 + s[0] - state_col.getJ(0)) < 2 &&
                                std::abs(l0 + s[0] - state_row.getJ(0)) < 2) {
                                array_j0.insert(l0 + s[0]);
                            }
                            for (float j0 : array_j0) {

                                std::set<float> array_j1;
                                if (std::abs(std::abs(l1 - s[1]) - state_col.getJ(1)) < 2 &&
                                    std::abs(std::abs(l1 - s[1]) - state_row.getJ(1)) < 2) {
                                    array_j1.insert(std::abs(l1 - s[1]));
                                }
                                if (std::abs(l1 + s[1] - state_col.getJ(1)) < 2 &&
                                    std::abs(l1 + s[1] - state_row.getJ(1)) < 2) {
                                    array_j1.insert(l1 + s[1]);
                                }
                                for (float j1 : array_j1) {

                                    for (auto &idx : array_q) {
                                        int q0_forth = idx[0]; // q = final.m - initial.m
                                        int q1_forth = idx[1];

                                        float m0 = state_col.getM(0) + q0_forth;
                                        if (std::abs(m0) > j0) {
                                            continue;
                                        }

                                        float m1 = state_col.getM(1) + q1_forth;
                                        if (std::abs(m1) > j1) {
                                            continue;
                                        }

                                        int q0_back = state_row.getM(0) - m0;
                                        int q1_back = state_row.getM(1) - m1;
                                        if (std::abs(q0_back) > 1 || std::abs(q1_back) > 1) {
                                            continue;
                                        }
                                        if (array_q.size() == 3 && q0_back + q1_back != 0) {
                                            continue;
                                        }

                                        StateTwo state_virtual =
                                            StateTwo(species, {{n0, n1}}, {{l0, l1}}, {{j0, j1}},
                                                     {{m0, m1}});
                                        if (set_states.find(state_virtual) != set_states.end()) {
                                            continue;
                                        }

                                        double energydiff_row =
                                            state_row.getEnergy() - state_virtual.getEnergy();
                                        if (bfield != 0) {
                                            energydiff_row += -bfield *
                                                (cache.getMagneticDipole(
                                                     state_row.getFirstState(),
                                                     state_row.getFirstState()) +
                                                 cache.getMagneticDipole(
                                                     state_row.getSecondState(),
                                                     state_row.getSecondState()) -
                                                 cache.getMagneticDipole(
                                                     state_virtual.getFirstState(),
                                                     state_virtual.getFirstState()) -
                                                 cache.getMagneticDipole(
                                                     state_virtual.getSecondState(),
                                                     state_virtual.getSecondState()));
                                        }

                                        double energydiff_col =
                                            state_col.getEnergy() - state_virtual.getEnergy();
                                        if (bfield != 0) {
                                            energydiff_col += -bfield *
                                                (cache.getMagneticDipole(
                                                     state_col.getFirstState(),
                                                     state_col.getFirstState()) +
                                                 cache.getMagneticDipole(
                                                     state_col.getSecondState(),
                                                     state_col.getSecondState()) -
                                                 cache.getMagneticDipole(
                                                     state_virtual.getFirstState(),
                                                     state_virtual.getFirstState()) -
                                                 cache.getMagneticDipole(
                                                     state_virtual.getSecondState(),
                                                     state_virtual.getSecondState()));
                                        }

                                        C6 += coulombs_constant *
                                            array_angle_term[3 * (q0_back + 1) + (q1_back + 1)] *
                                            cache.getElectricDipole(state_row.getFirstState(),
                                                                    state_virtual.getFirstState()) *
                                            cache.getElectricDipole(
                                                state_row.getSecondState(),
                                                state_virtual.getSecondState()) *
                                            coulombs_constant *
                                            array_angle_term[3 * (q0_forth + 1) + (q1_forth + 1)] *
                                            cache.getElectricDipole(state_virtual.getFirstState(),
                                                                    state_col.getFirstState()) *
                                            cache.getElectricDipole(state_virtual.getSecondState(),
                                                                    state_col.getSecondState()) *
                                            0.5 *
                                            (1 / energydiff_row +
                                             1 / energydiff_col); // getMultipole(final, inital, 1)
                                    }
                                }
                            }
                        }
                    }
                }
            }

            C6_matrix(idx_row, idx_col) = C6;
        }
    }

    return C6_matrix.selfadjointView<Eigen::Upper>();
}

eigen_dense_double_t PerturbativeInteraction::getC3(const std::vector<StateTwo> &states) {
    eigen_dense_double_t C3_matrix = eigen_dense_double_t::Zero(states.size(), states.size());

    for (size_t idx_row = 0; idx_row < states.size(); ++idx_row) {
        auto &state_row = states[idx_row];

        for (size_t idx_col = idx_row + 1; idx_col < states.size(); ++idx_col) {
            auto &state_col = states[idx_col];

            int q0 = state_row.getM(0) - state_col.getM(0);
            int q1 = state_row.getM(1) - state_col.getM(1);

            if (array_q.size() == 3 && q0 + q1 != 0) {
                continue;
            }
            if (std::abs(q0) > 1 || std::abs(q1) > 1) {
                continue;
            }

            C3_matrix(idx_row, idx_col) = coulombs_constant *
                array_angle_term[3 * (q0 + 1) + (q1 + 1)] *
                cache.getElectricDipole(state_row.getFirstState(), state_col.getFirstState()) *
                cache.getElectricDipole(state_row.getSecondState(), state_col.getSecondState());
        }
    }

    return C3_matrix.selfadjointView<Eigen::Upper>();
}

eigen_dense_double_t PerturbativeInteraction::getEnergy(const std::vector<StateTwo> &states) {
    eigen_dense_double_t energies_matrix = eigen_dense_double_t::Zero(states.size(), states.size());

    for (size_t idx = 0; idx < states.size(); ++idx) {
        auto &state = states[idx];

        double energy = state.getEnergy();
        if (bfield != 0) {
            energy += -bfield *
                (cache.getMagneticDipole(state.getFirstState(), state.getFirstState()) +
                 cache.getMagneticDipole(state.getSecondState(), state.getSecondState()));
        }
        energies_matrix(idx, idx) = energy;
    }

    return energies_matrix;
}

void PerturbativeInteraction::initializeAngleTerms(double angle) {
    if (angle == 0) {
        array_q.reserve(3);

        array_q.push_back({{0, 0}});
        array_angle_term[3 * (0 + 1) + (0 + 1)] = -2;

        array_q.push_back({{1, -1}});
        array_angle_term[3 * (1 + 1) + (-1 + 1)] = -1;

        array_q.push_back({{-1, 1}});
        array_angle_term[3 * (-1 + 1) + (1 + 1)] = -1;
    } else {
        array_q.reserve(9);

        array_q.push_back({{0, 0}});
        array_angle_term[3 * (0 + 1) + (0 + 1)] = 1. - 3. * std::pow(std::cos(angle), 2);

        array_q.push_back({{1, -1}});
        array_angle_term[3 * (1 + 1) + (-1 + 1)] = -1 + 1.5 * std::pow(std::sin(angle), 2);

        array_q.push_back({{-1, 1}});
        array_angle_term[3 * (-1 + 1) + (1 + 1)] = -1 + 1.5 * std::pow(std::sin(angle), 2);

        array_q.push_back({{1, 1}});
        array_angle_term[3 * (1 + 1) + (1 + 1)] = -1.5 * std::pow(std::sin(angle), 2);

        array_q.push_back({{-1, -1}});
        array_angle_term[3 * (-1 + 1) + (-1 + 1)] = -1.5 * std::pow(std::sin(angle), 2);

        array_q.push_back({{0, 1}});
        array_angle_term[3 * (0 + 1) + (1 + 1)] =
            3. / std::sqrt(2) * std::sin(angle) * std::cos(angle);

        array_q.push_back({{1, 0}});
        array_angle_term[3 * (1 + 1) + (0 + 1)] =
            3. / std::sqrt(2) * std::sin(angle) * std::cos(angle);

        array_q.push_back({{0, -1}});
        array_angle_term[3 * (0 + 1) + (-1 + 1)] =
            -3. / std::sqrt(2) * std::sin(angle) * std::cos(angle);

        array_q.push_back({{-1, 0}});
        array_angle_term[3 * (-1 + 1) + (0 + 1)] =
            -3. / std::sqrt(2) * std::sin(angle) * std::cos(angle);
    }
}
