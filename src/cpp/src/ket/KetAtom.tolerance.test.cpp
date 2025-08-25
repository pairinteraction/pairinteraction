// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/Parity.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {

DOCTEST_TEST_CASE("create ket with small numerical errors in quantum numbers") {
    Database &database = Database::get_global_instance();
    
    // Test case: Create a KetAtom with quantum numbers that have small numerical errors
    // This simulates what happens when values from database are like 2.9999999 instead of 3.0
    
    double quantum_number_j_exp_with_error = 3.0 - 1e-12;  // Very close to 3 but not exactly
    double quantum_number_l_exp_with_error = 3.0 - 1e-12;  // Very close to 3 but not exactly
    double quantum_number_f = 3.0 - 1e-12;  // Very close to 3 but not exactly
    double quantum_number_m = 0.5;
    
    // This should not fail due to the numerical tolerance fixes
    DOCTEST_CHECK_NOTHROW({
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
        
        // Try to get the label - this is where the abort() would have happened before the fix
        std::string label = ket.get_label();
        DOCTEST_CHECK(!label.empty());
        DOCTEST_MESSAGE("Generated label with numerical tolerance: ", label);
    });
}

DOCTEST_TEST_CASE("create ket with larger numerical errors should still work") {
    Database &database = Database::get_global_instance();
    
    // Test with slightly larger but still reasonable numerical errors
    double quantum_number_j_exp_with_error = 2.0 + 1e-10;  // Very close to 2 but not exactly
    double quantum_number_l_exp_with_error = 1.0 - 1e-10;  // Very close to 1 but not exactly
    double quantum_number_f = 2.0 + 1e-10;  // Very close to 2 but not exactly
    double quantum_number_m = 0.5 - 1e-12;  // Very close to 0.5 but not exactly
    
    DOCTEST_CHECK_NOTHROW({
        KetAtom ket(typename KetAtom::Private(), 
                   -1.0,  // energy
                   quantum_number_f,  // f
                   quantum_number_m,  // m
                   Parity::ODD,  // parity
                   "Rb",  // species
                   60,  // n
                   60.0,  // nu
                   0.0,  // nui_exp
                   0.0,  // nui_std
                   quantum_number_l_exp_with_error,  // l_exp - this has numerical error
                   0.0,  // l_std
                   0.5,  // s_exp
                   0.0,  // s_std
                   quantum_number_j_exp_with_error,  // j_exp - this has numerical error
                   0.0,  // j_std
                   0.0,  // l_ryd_exp
                   0.0,  // l_ryd_std
                   0.0,  // j_ryd_exp
                   0.0,  // j_ryd_std
                   false,  // is_j_total_momentum
                   false,  // is_calculated_with_mqdt
                   0.0,  // underspecified_channel_contribution
                   database,
                   456);  // id_in_database
        
        std::string label = ket.get_label();
        DOCTEST_CHECK(!label.empty());
        DOCTEST_MESSAGE("Generated label with larger numerical tolerance: ", label);
    });
}

} // namespace pairinteraction