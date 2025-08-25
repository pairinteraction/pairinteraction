// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

// Test to reproduce the numerical precision issue in KetAtom
// This test demonstrates the problem where KetAtom fails with std::abort()
// when quantum numbers have small numerical errors (e.g., 2.9999999 instead of 3.0)

#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/Parity.hpp"

#include <iostream>
#include <limits>

namespace pairinteraction {

void test_numerical_precision_issue() {
    std::cout << "Testing numerical precision issue in KetAtom...\n";
    
    // Create a mock database - this is just for testing the KetAtom constructor directly
    Database database;
    
    // Test case: Create a KetAtom with quantum numbers that have small numerical errors
    // This simulates what happens when values from database are like 2.9999999 instead of 3.0
    
    double quantum_number_j_exp_with_error = 3.0 - 1e-12;  // Very close to 3 but not exactly
    double quantum_number_l_exp_with_error = 3.0 - 1e-12;  // Very close to 3 but not exactly
    double quantum_number_f = 3.0 - 1e-12;  // Very close to 3 but not exactly
    double quantum_number_m = 0.5;
    
    try {
        // This should fail before the fix due to exact equality checks in get_label()
        KetAtom ket(typename KetAtom::Private(), 
                   -1.0,  // energy
                   quantum_number_f,  // f
                   quantum_number_m,  // m
                   Parity::ODD,  // parity
                   "Yb174_mqdt",  // species
                   30,  // n
                   30.0,  // nu
                   0.0,  // nui_exp
                   0.0,  // nui_std
                   quantum_number_l_exp_with_error,  // l_exp - this has numerical error
                   0.0,  // l_std
                   0.5,  // s_exp
                   0.0,  // s_std
                   quantum_number_j_exp_with_error,  // j_exp - this has numerical error
                   0.0,  // j_std
                   quantum_number_l_exp_with_error,  // l_ryd_exp - this has numerical error
                   0.0,  // l_ryd_std
                   quantum_number_j_exp_with_error,  // j_ryd_exp - this has numerical error
                   0.0,  // j_ryd_std
                   true,  // is_j_total_momentum
                   true,  // is_calculated_with_mqdt
                   0.0,  // underspecified_channel_contribution
                   database,
                   123);  // id_in_database
        
        // Try to get the label - this is where the abort() would happen
        std::string label = ket.get_label();
        std::cout << "Success! Generated label: " << label << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
    }
}

} // namespace pairinteraction

int main() {
    pairinteraction::test_numerical_precision_issue();
    return 0;
}