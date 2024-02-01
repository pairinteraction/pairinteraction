#include <streambuf>
#include <sstream>

#include "Constants.hpp"
#include "Interface.hpp"
#include "MatrixElementCache.hpp"
#include "PerturbativeInteraction.hpp"
#include "State.hpp"
#include "Symmetry.hpp"
#include "SystemBase.hpp"
#include "SystemOne.hpp"
#include "SystemTwo.hpp"
#include "QuantumDefect.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/eigen/dense.h>
//#include <nanobind/eigen/sparse.h>

namespace nb = nanobind;

using namespace nb::literals;

struct array_source : std::streambuf {
    array_source(char *buffer, size_t len) {
        this->setg(buffer, buffer, buffer + len);
    }
};

template <typename T>
PyObject* __getstate__(T const &self, std::string const &cls) {
    std::stringstream ss;
    cereal::BinaryOutputArchive ar(ss);
    ar << cereal::make_nvp(cls, self);
    return PyBytes_FromStringAndSize(ss.str().data(), ss.str().length());
}

template <typename T>
void __setstate__(PyObject* const sState, T &self, std::string const &cls) {
    char *buffer;
    Py_ssize_t len;
    PyBytes_AsStringAndSize(sState, &buffer, &len);
    array_source asource(buffer, static_cast<size_t>(len));
    std::istream ss(&asource);
    cereal::BinaryInputArchive ar(ss);
    ar >> cereal::make_nvp(cls, self);
}

template <typename System, typename NBClass>
void add_system_base_methods(NBClass &system) {
    using State = typename System::State;
    system
        .def("setMinimalNorm", &System::setMinimalNorm)
        .def("restrictEnergy", &System::restrictEnergy)

        .def("restrictN", nb::overload_cast<int, int>(&System::restrictN))
        .def("restrictN", nb::overload_cast<std::set<int>>(&System::restrictN))
        .def("restrictL", nb::overload_cast<int, int>(&System::restrictL))
        .def("restrictL", nb::overload_cast<std::set<int>>(&System::restrictL))
        .def("restrictJ", nb::overload_cast<float, float>(&System::restrictJ))
        .def("restrictJ", nb::overload_cast<std::set<float>>(&System::restrictJ))
        .def("restrictM", nb::overload_cast<float, float>(&System::restrictM))
        .def("restrictM", nb::overload_cast<std::set<float>>(&System::restrictM))

        .def("addStates", nb::overload_cast<const State &>(&System::addStates))
        .def("addStates", nb::overload_cast<const std::set<State> &>(&System::addStates))

        .def("getOverlap", nb::overload_cast<const State &>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const std::vector<State> &>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const size_t &>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const std::vector<size_t> &>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const State &, std::array<double, 3>, std::array<double, 3>>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const std::vector<State> &, std::array<double, 3>, std::array<double, 3>>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const size_t &, std::array<double, 3>, std::array<double, 3>>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const std::vector<size_t> &, std::array<double, 3>, std::array<double, 3>>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const State &, double, double, double>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const size_t &, double, double, double>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const std::vector<State> &, double, double, double>(&System::getOverlap))
        .def("getOverlap", nb::overload_cast<const std::vector<size_t> &, double, double, double>(&System::getOverlap))

        .def("getStates", &System::getStates)
        .def("getStatesMultiIndex", &System::getStatesMultiIndex)
        .def("getBasisvectors", &System::getBasisvectors)
        .def("getHamiltonian",  &System::getHamiltonian)
        .def("getNumBasisvectors", &System::getNumBasisvectors)
        .def("getNumStates", &System::getNumStates)
        .def("getMainStates", &System::getMainStates)
        .def("getConnections", &System::getConnections)

        .def("buildHamiltonian", &System::buildHamiltonian)
        .def("buildInteraction", &System::buildInteraction)
        .def("buildBasis", &System::buildBasis)
        .def("diagonalize", nb::overload_cast<double, double>(&System::diagonalize))
        .def("diagonalize", nb::overload_cast<double, double, double>(&System::diagonalize))
        .def("diagonalize", nb::overload_cast<>(&System::diagonalize))
        .def("diagonalize", nb::overload_cast<double>(&System::diagonalize))
        .def("canonicalize", &System::canonicalize)
        .def("unitarize", &System::unitarize)
        .def("rotate", nb::overload_cast<std::array<double, 3>, std::array<double, 3>>(&System::rotate))
        .def("rotate", nb::overload_cast<double, double, double>(&System::rotate))
        .def("add", &System::add)
        .def("constrainBasisvectors", &System::constrainBasisvectors)
        .def("applySchriefferWolffTransformation", &System::applySchriefferWolffTransformation)

        .def("getStateIndex", nb::overload_cast<const State &>(&System::getStateIndex))
        .def("getStateIndex", nb::overload_cast<const std::vector<State> &>(&System::getStateIndex))
        .def("getBasisvectorIndex", nb::overload_cast<const State &>(&System::getBasisvectorIndex))
        .def("getBasisvectorIndex", nb::overload_cast<const std::vector<State> &>(&System::getBasisvectorIndex))
        .def("forgetStatemixing", &System::forgetStatemixing)
        .def("getHamiltonianEntry", &System::getHamiltonianEntry)
        .def("setHamiltonianEntry", &System::setHamiltonianEntry)
        .def("addHamiltonianEntry", &System::addHamiltonianEntry);
}

template <typename Scalar>
void declare_systems(nb::module_ &m, std::string const &type) {
    std::string pyclass_name_system_one_base = "_SystemStateOne" + type;
    nb::class_<SystemBase<Scalar, StateOne>> system_one_base(m, pyclass_name_system_one_base.c_str());
    add_system_base_methods<SystemBase<Scalar, StateOne>>(system_one_base);

    std::string pyclass_name_system_one = "SystemOne" + type;
    nb::class_<SystemOne<Scalar>, SystemBase<Scalar, StateOne>> system_one(m, pyclass_name_system_one.c_str());
    add_system_base_methods<SystemOne<Scalar>>(system_one);
    system_one
        .def(nb::init<SystemOne<Scalar> const &>())
        .def(nb::init<std::string, MatrixElementCache &>())
        .def(nb::init<std::string, MatrixElementCache &, bool>())
        .def("getSpecies", &SystemOne<Scalar>::getSpecies)
        .def("setEfield", nb::overload_cast<std::array<double, 3>>(&SystemOne<Scalar>::setEfield))
        .def("setEfield", nb::overload_cast<std::array<double, 3>, std::array<double, 3>, std::array<double, 3>>(&SystemOne<Scalar>::setEfield))
        .def("setEfield", nb::overload_cast<std::array<double, 3>, double, double, double>(&SystemOne<Scalar>::setEfield))
        .def("setBfield", nb::overload_cast<std::array<double, 3>>(&SystemOne<Scalar>::setBfield))
        .def("setBfield", nb::overload_cast<std::array<double, 3>, std::array<double, 3>, std::array<double, 3>>(&SystemOne<Scalar>::setBfield))
        .def("setBfield", nb::overload_cast<std::array<double, 3>, double, double, double>(&SystemOne<Scalar>::setBfield))
        .def("enableDiamagnetism", &SystemOne<Scalar>::enableDiamagnetism)
        .def("setIonCharge", &SystemOne<Scalar>::setIonCharge)
        .def("setRydIonOrder", &SystemOne<Scalar>::setRydIonOrder)
        .def("setRydIonDistance", &SystemOne<Scalar>::setRydIonDistance)
        .def("setConservedParityUnderReflection", &SystemOne<Scalar>::setConservedParityUnderReflection)
        .def("setConservedMomentaUnderRotation", &SystemOne<Scalar>::setConservedMomentaUnderRotation)
        .def("__getstate__", [pyclass_name_system_one] (SystemOne<Scalar> const &self) -> nb::handle {
            return __getstate__(self, pyclass_name_system_one);
        })
        .def("__setstate__", [pyclass_name_system_one] (SystemOne<Scalar> &self, nb::handle buf) {
            MatrixElementCache tmp;
            new (&self) SystemOne<Scalar>("", tmp);
            __setstate__(buf.ptr(), self, pyclass_name_system_one);
        });

    std::string pyclass_name_system_two_base = "_SystemStateTwo" + type;
    nb::class_<SystemBase<Scalar, StateTwo>> system_two_base(m, pyclass_name_system_two_base.c_str());
    add_system_base_methods<SystemBase<Scalar, StateTwo>>(system_two_base);

    std::string pyclass_name_system_two = "SystemTwo" + type;
    nb::class_<SystemTwo<Scalar>, SystemBase<Scalar, StateTwo>> system_two(m, pyclass_name_system_two.c_str());
    add_system_base_methods<SystemTwo<Scalar>>(system_two);
    system_two
        .def(nb::init<SystemTwo<Scalar> const &>())
        .def(nb::init<const SystemOne<Scalar> &, const SystemOne<Scalar> &, MatrixElementCache &>())
        .def(nb::init<const SystemOne<Scalar> &, const SystemOne<Scalar> &, MatrixElementCache &, bool>())
        .def("getSpecies", &SystemTwo<Scalar>::getSpecies)
        .def("getStatesFirst", &SystemTwo<Scalar>::getStatesFirst)
        .def("getStatesSecond", &SystemTwo<Scalar>::getStatesSecond)
        .def("enableGreenTensor", &SystemTwo<Scalar>::enableGreenTensor)
        .def("setSurfaceDistance", &SystemTwo<Scalar>::setSurfaceDistance)
        .def("setAngle", &SystemTwo<Scalar>::setAngle)
        .def("setDistance", &SystemTwo<Scalar>::setDistance)
        .def("setDistanceVector", &SystemTwo<Scalar>::setDistanceVector)
        .def("setOrder", &SystemTwo<Scalar>::setOrder)
        .def("setConservedParityUnderPermutation", &SystemTwo<Scalar>::setConservedParityUnderPermutation)
        .def("setConservedParityUnderInversion", &SystemTwo<Scalar>::setConservedParityUnderInversion)
        .def("setConservedParityUnderReflection", &SystemTwo<Scalar>::setConservedParityUnderReflection)
        .def("setConservedMomentaUnderRotation", &SystemTwo<Scalar>::setConservedMomentaUnderRotation)
        .def("setOneAtomBasisvectors", &SystemTwo<Scalar>::setOneAtomBasisvectors)
        .def("__getstate__", [pyclass_name_system_two] (SystemTwo<Scalar> const &self) -> nb::handle {
            return __getstate__(self, pyclass_name_system_two);
        })
        .def("__setstate__", [pyclass_name_system_two] (SystemTwo<Scalar> &self, nb::handle buf) {
            MatrixElementCache tmp;
            SystemOne<Scalar> s1("", tmp);
            SystemOne<Scalar> s2("", tmp);
            new (&self) SystemTwo<Scalar>(s1, s2, tmp);
            __setstate__(buf.ptr(), self, pyclass_name_system_two);
        });
}

NB_MODULE(binding, m) {

    // Constants.hpp

    m.attr("au2GHz") = au2GHz;
    m.attr("au2Vcm") = au2Vcm;
    m.attr("au2G") = au2G;
    m.attr("au2um") = au2um;
    m.attr("inverse_electric_constant") = inverse_electric_constant;
    m.attr("sqrt_inverse_electric_constant") = sqrt_inverse_electric_constant;
    m.attr("inverse_electron_rest_mass") = inverse_electron_rest_mass;
    m.attr("coulombs_constant") = coulombs_constant;
    m.attr("electron_rest_mass") = electron_rest_mass;
    m.attr("elementary_charge") = elementary_charge;
    m.attr("bohr_magneton") = bohr_magneton;
    m.attr("reduced_planck_constant") = reduced_planck_constant;
    m.attr("speed_of_light") = speed_of_light;
    m.attr("muB") = muB;
    m.attr("gS") = gS;
    m.attr("gL") = gL;
    m.attr("ARB") = ARB;
    m.attr("mkl_enabled") = mkl_enabled;
    m.attr("gsl_enabled") = gsl_enabled;

    // Interface.hpp

    m.def("computeReal", &compute<double>);
    m.def("computeComplex", &compute<std::complex<double>>);

    // Symmetry.hpp

    nb::enum_<parity_t>(m, "parity")
        .value("NA", NA)
        .value("EVEN", EVEN)
        .value("ODD", ODD)
        .export_values();

    nb::class_<Symmetry> symmetry(m, "Symmetry");
    symmetry
        .def(nb::init<>())
        .def_rw("inversion", &Symmetry::inversion)
        .def_rw("reflection", &Symmetry::reflection)
        .def_rw("permutation", &Symmetry::permutation)
        .def_rw("rotation", &Symmetry::rotation)
        .def(nb::self < nb::self);

    // MatrixElementCache.hpp

    nb::enum_<method_t>(m, "method")
        .value("NUMEROV", NUMEROV)
        .value("WHITTAKER", WHITTAKER)
        .export_values();

    m.def("selectionRulesMomentumNew", nb::overload_cast<StateOne const &, StateOne const &, int>(&selectionRulesMomentumNew));
    m.def("selectionRulesMomentumNew", nb::overload_cast<StateOne const &, StateOne const &>(&selectionRulesMomentumNew));
    m.def("selectionRulesMultipoleNew", nb::overload_cast<StateOne const &, StateOne const &, int, int>(&selectionRulesMultipoleNew));
    m.def("selectionRulesMultipoleNew", nb::overload_cast<StateOne const &, StateOne const &, int>(&selectionRulesMultipoleNew));

    nb::class_<MatrixElementCache> matrix_element_cache(m, "MatrixElementCache");
    matrix_element_cache
        .def(nb::init<>())
        .def(nb::init<std::string>())
        .def("getElectricDipole", &MatrixElementCache::getElectricDipole)
        .def("getElectricMultipole", nb::overload_cast<const StateOne &, const StateOne &, int>(&MatrixElementCache::getElectricMultipole))
        .def("getDiamagnetism", &MatrixElementCache::getDiamagnetism)
        .def("getMagneticDipole", &MatrixElementCache::getMagneticDipole)
        .def("getElectricMultipole", nb::overload_cast<const StateOne &, const StateOne &, int, int>(&MatrixElementCache::getElectricMultipole))
        .def("getRadial", &MatrixElementCache::getRadial)
        .def("precalculateElectricMomentum", &MatrixElementCache::precalculateElectricMomentum)
        .def("precalculateMagneticMomentum", &MatrixElementCache::precalculateMagneticMomentum)
        .def("precalculateDiamagnetism", &MatrixElementCache::precalculateDiamagnetism)
        .def("precalculateMultipole", &MatrixElementCache::precalculateMultipole)
        .def("precalculateRadial", &MatrixElementCache::precalculateRadial)
        .def("setDefectDB", &MatrixElementCache::setDefectDB)
        .def("getDefectDB", &MatrixElementCache::getDefectDB)
        .def("setMethod", &MatrixElementCache::setMethod)
        .def("loadElectricDipoleDB", &MatrixElementCache::loadElectricDipoleDB)
        .def("size", &MatrixElementCache::size)
        .def("__getstate__", [] (MatrixElementCache const &self) -> nb::handle {
            return __getstate__(self, "MatrixElementCache");
        })
        .def("__setstate__", [] (MatrixElementCache &self, nb::handle buf) {
            new (&self) MatrixElementCache;
            __setstate__(buf.ptr(), self, "MatrixElementCache");
        });

    // State.hpp

    nb::class_<StateOne> state_one(m, "StateOne");
    state_one
        .def(nb::init<std::string, int, int, float, float>())
        .def(nb::init<std::string>())
        .def("__str__", &StateOne::str)
        .def("getN", &StateOne::getN)
        .def("getL", &StateOne::getL)
        .def("getJ", &StateOne::getJ)
        .def("getM", &StateOne::getM)
        .def("getS", &StateOne::getS)
        .def("getSpecies", &StateOne::getSpecies)
        .def("getElement", &StateOne::getElement)
        .def("getEnergy", nb::overload_cast<>(&StateOne::getEnergy, nb::const_))
        .def("getEnergy", nb::overload_cast<MatrixElementCache &>(&StateOne::getEnergy, nb::const_))
        .def("getNStar", nb::overload_cast<>(&StateOne::getNStar, nb::const_))
        .def("getNStar", nb::overload_cast<MatrixElementCache &>(&StateOne::getNStar, nb::const_))
        .def("getLabel", &StateOne::getLabel)
        .def("isArtificial", &StateOne::isArtificial)
        .def("isGeneralized", &StateOne::isGeneralized)
        .def("getHash", &StateOne::getHash)
        .def("getReflected", &StateOne::getReflected)
        .def(nb::self == nb::self)
        .def(nb::self ^ nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        .def(nb::self <= nb::self)
        .def("__getstate__", [] (StateOne const &self) -> nb::handle {
            return __getstate__(self, "StateOne");
        })
        .def("__setstate__", [] (StateOne &self, nb::handle buf) {
            new (&self) StateOne;
            __setstate__(buf.ptr(), self, "StateOne");
        });


    nb::class_<StateTwo> state_two(m, "StateTwo");
    state_two
        .def(nb::init<std::array<std::string, 2>, std::array<int, 2>,
             std::array<int, 2>, std::array<float, 2>, std::array<float, 2>>())
        .def(nb::init<std::array<std::string, 2>>())
        .def(nb::init<StateOne, StateOne>())
        .def("__str__", &StateTwo::str)
        .def("getN", nb::overload_cast<>(&StateTwo::getN, nb::const_))
        .def("getL", nb::overload_cast<>(&StateTwo::getL, nb::const_))
        .def("getJ", nb::overload_cast<>(&StateTwo::getJ, nb::const_))
        .def("getM", nb::overload_cast<>(&StateTwo::getM, nb::const_))
        .def("getS", nb::overload_cast<>(&StateTwo::getS, nb::const_))
        .def("getSpecies", nb::overload_cast<>(&StateTwo::getSpecies, nb::const_))
        .def("getElement", nb::overload_cast<>(&StateTwo::getElement, nb::const_))
        .def("getEnergy", nb::overload_cast<>(&StateTwo::getEnergy, nb::const_))
        .def("getEnergy", nb::overload_cast<MatrixElementCache &>(&StateTwo::getEnergy, nb::const_))
        .def("getNStar", nb::overload_cast<>(&StateTwo::getNStar, nb::const_))
        .def("getNStar", nb::overload_cast<MatrixElementCache &>(&StateTwo::getNStar, nb::const_))
        .def("getLeRoyRadius", &StateTwo::getLeRoyRadius)
        .def("getLabel", nb::overload_cast<>(&StateTwo::getLabel, nb::const_))
        .def("isArtificial", nb::overload_cast<>(&StateTwo::isArtificial, nb::const_))
        .def("isGeneralized", nb::overload_cast<>(&StateTwo::isGeneralized, nb::const_))
        .def("getN", nb::overload_cast<int>(&StateTwo::getN, nb::const_))
        .def("getL", nb::overload_cast<int>(&StateTwo::getL, nb::const_))
        .def("getJ", nb::overload_cast<int>(&StateTwo::getJ, nb::const_))
        .def("getM", nb::overload_cast<int>(&StateTwo::getM, nb::const_))
        .def("getS", nb::overload_cast<int>(&StateTwo::getS, nb::const_))
        .def("getSpecies", nb::overload_cast<int>(&StateTwo::getSpecies, nb::const_))
        .def("getElement", nb::overload_cast<int>(&StateTwo::getElement, nb::const_))
        .def("getEnergy", nb::overload_cast<int>(&StateTwo::getEnergy, nb::const_))
        .def("getEnergy", nb::overload_cast<int, MatrixElementCache &>(&StateTwo::getEnergy, nb::const_))
        .def("getNStar", nb::overload_cast<int>(&StateTwo::getNStar, nb::const_))
        .def("getNStar", nb::overload_cast<int, MatrixElementCache &>(&StateTwo::getNStar, nb::const_))
        .def("getLabel", nb::overload_cast<int>(&StateTwo::getLabel, nb::const_))
        .def("isArtificial", nb::overload_cast<int>(&StateTwo::isArtificial, nb::const_))
        .def("isGeneralized", nb::overload_cast<int>(&StateTwo::isGeneralized, nb::const_))
        .def("getFirstState", &StateTwo::getFirstState)
        .def("getSecondState", &StateTwo::getSecondState)
        .def("getHash", &StateTwo::getHash)
        .def("getReflected", &StateTwo::getReflected)
        .def(nb::self == nb::self)
        .def(nb::self ^ nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        .def(nb::self <= nb::self)
        .def("__getstate__", [] (StateTwo const &self) -> nb::handle {
            return __getstate__(self, "StateTwo");
        })
        .def("__setstate__", [] (StateTwo &self, nb::handle buf) {
            new (&self) StateTwo;
            __setstate__(buf.ptr(), self, "StateTwo");
        });


    // QuantumDefect.hpp

    nb::class_<QuantumDefect> quantum_defect(m, "QuantumDefect");
    quantum_defect
        .def(nb::init<std::string const &, int, int, double>())
        .def(nb::init<std::string const &, int, int, double, std::string const &>())
        .def_prop_ro("species", [](QuantumDefect &qd) { return qd.species; })
        .def_prop_ro("n", [](QuantumDefect &qd) { return qd.n; })
        .def_prop_ro("l", [](QuantumDefect &qd) { return qd.l; })
        .def_prop_ro("j", [](QuantumDefect &qd) { return qd.j; })
        .def_prop_ro("ac", [](QuantumDefect &qd) { return qd.ac; })
        .def_prop_ro("Z", [](QuantumDefect &qd) { return qd.Z; })
        .def_prop_ro("a1", [](QuantumDefect &qd) { return qd.a1; })
        .def_prop_ro("a2", [](QuantumDefect &qd) { return qd.a2; })
        .def_prop_ro("a3", [](QuantumDefect &qd) { return qd.a3; })
        .def_prop_ro("a4", [](QuantumDefect &qd) { return qd.a4; })
        .def_prop_ro("rc", [](QuantumDefect &qd) { return qd.rc; })
        .def_prop_ro("nstar", [](QuantumDefect &qd) { return qd.nstar; })
        .def_prop_ro("energy", [](QuantumDefect &qd) { return qd.energy; });

    m.def("energy_level", &energy_level);
    m.def("nstar", &nstar);

    // Wavefunction.hpp

    m.def("V", &model_potential::V);
    m.def("g", &model_potential::g);

    nb::class_<Numerov> numerov(m, "Numerov");
    numerov
        .def(nb::init<QuantumDefect const &>())
        .def_ro_static("dx", &Numerov::dx)
        .def("integrate", &Numerov::integrate)
        .def_static("power_kernel", &Numerov::power_kernel);

    m.def("HypergeometricU", &whittaker_functions::HypergeometricU);
    m.def("WhittakerW", &whittaker_functions::WhittakerW);
    m.def("RadialWFWhittaker", &whittaker_functions::RadialWFWhittaker);

    nb::class_<Whittaker> whittaker(m, "Whittaker");
    whittaker
        .def(nb::init<QuantumDefect const &>())
        .def_ro_static("dx", &Whittaker::dx)
        .def("integrate", &Whittaker::integrate)
        .def_static("power_kernel", &Whittaker::power_kernel);

    // PerturbativeInteraction.hpp

    nb::class_<PerturbativeInteraction> perturbative_interaction(m, "PerturbativeInteraction");
    perturbative_interaction
        .def(nb::init<MatrixElementCache &>())
        .def(nb::init<double, MatrixElementCache &>())
        .def("getC6", nb::overload_cast<const StateTwo &, double>(&PerturbativeInteraction::getC6))
        .def("getC6", nb::overload_cast<const std::vector<StateTwo> &, double>(&PerturbativeInteraction::getC6))
        .def("getC3", &PerturbativeInteraction::getC3)
        .def("getEnergy", &PerturbativeInteraction::getEnergy);

    // WignerD.hpp

    nb::class_<WignerD> wigner_d(m, "WignerD");
    wigner_d
        .def(nb::init<>())
        .def("__call__", nb::overload_cast<float, float, float, double>(&WignerD::operator()))
        .def("__call__", nb::overload_cast<float, float, float, double, double, double>(&WignerD::operator()));

    // SystemOne.hpp SystemTwo.hpp
    declare_systems<double>(m, "Real");
    declare_systems<std::complex<double>>(m, "Complex");
}
