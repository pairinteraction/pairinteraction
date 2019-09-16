#include <Eigen/Sparse>

#include <jlcxx/jlcxx.hpp>
#include <jlcxx/const_array.hpp>
#include "State.h"
#include "SystemBase.h"
#include "SystemOne.h"
#include "SystemTwo.h"
#include "MatrixElementCache.h"
#include "dtypes.h"

namespace jlcxx
{
  template<> struct IsBits<method_t> : std::true_type {};
  template<> struct IsBits<parity_t> : std::true_type {};
}

JLCXX_MODULE define_julia_module(jlcxx::Module& pi)
{
  pi.add_bits<method_t>("method_t");
  pi.set_const("NUMEROV", NUMEROV);
  pi.set_const("WHITTAKER", WHITTAKER);

  pi.add_bits<parity_t>("parity_t");
  pi.set_const("NA", NA);
  pi.set_const("EVEN", EVEN);
  pi.set_const("ODD",   ODD);

  pi.set_const("ARB", ARB);

  pi.add_type<MatrixElementCache>("MatrixElementCache")
    .constructor<std::string>()
    .method("getElectricDipole", &MatrixElementCache::getElectricDipole)
    .method("getElectricMultipole", static_cast<double (MatrixElementCache::*)(const StateOne&, const StateOne&, int)>(&MatrixElementCache::getElectricMultipole))
    .method("getDiamagnetism", &MatrixElementCache::getDiamagnetism)
    .method("getMagneticDipole", &MatrixElementCache::getMagneticDipole)
    .method("getElectricMultipole", static_cast<double (MatrixElementCache::*)(const StateOne&, const StateOne&, int, int)>(&MatrixElementCache::getElectricMultipole))
    .method("getRadial", &MatrixElementCache::getRadial)
    .method("getLeRoyRadius", &MatrixElementCache::getLeRoyRadius)
    .method("precalculateElectricMomentum", [](MatrixElementCache& mec, jlcxx::ArrayRef<jl_value_t*> basis_one_jl, int q) {
      std::vector<StateOne> basis_one;
      for (unsigned i=0; i<basis_one_jl.size(); i++){
        const StateOne s = *jlcxx::unbox_wrapped_ptr<StateOne>(basis_one_jl[i]);
        basis_one.push_back(s);
      }
      const std::vector<StateOne>& basis_one_const = basis_one;
      mec.precalculateElectricMomentum(basis_one_const, q);
    })
    .method("precalculateMagneticMomentum", [](MatrixElementCache& mec, jlcxx::ArrayRef<jl_value_t*> basis_one_jl, int q) {
      std::vector<StateOne> basis_one;
      for (unsigned i=0; i<basis_one_jl.size(); i++){
        const StateOne s = *jlcxx::unbox_wrapped_ptr<StateOne>(basis_one_jl[i]);
        basis_one.push_back(s);
      }
      const std::vector<StateOne>& basis_one_const = basis_one;
      mec.precalculateMagneticMomentum(basis_one_const, q);
    })
    .method("precalculateDiamagnetism", [](MatrixElementCache& mec, jlcxx::ArrayRef<jl_value_t*> basis_one_jl, int k, int q) {
      std::vector<StateOne> basis_one;
      for (unsigned i=0; i<basis_one_jl.size(); i++){
        const StateOne s = *jlcxx::unbox_wrapped_ptr<StateOne>(basis_one_jl[i]);
        basis_one.push_back(s);
      }
      const std::vector<StateOne>& basis_one_const = basis_one;
      mec.precalculateDiamagnetism(basis_one_const, k, q);
    })
    .method("precalculateMultipole", [](MatrixElementCache& mec, jlcxx::ArrayRef<jl_value_t*> basis_one_jl, int k) {
      std::vector<StateOne> basis_one;
      for (unsigned i=0; i<basis_one_jl.size(); i++){
        const StateOne s = *jlcxx::unbox_wrapped_ptr<StateOne>(basis_one_jl[i]);
        basis_one.push_back(s);
      }
      const std::vector<StateOne>& basis_one_const = basis_one;
      mec.precalculateMultipole(basis_one_const, k);
    })
    .method("precalculateRadial", [](MatrixElementCache& mec, jlcxx::ArrayRef<jl_value_t*> basis_one_jl, int k) {
      std::vector<StateOne> basis_one;
      for (unsigned i=0; i<basis_one_jl.size(); i++){
        const StateOne s = *jlcxx::unbox_wrapped_ptr<StateOne>(basis_one_jl[i]);
        basis_one.push_back(s);
      }
      const std::vector<StateOne>& basis_one_const = basis_one;
      mec.precalculateRadial(basis_one_const, k);
    })
    .method("setDefectDB", &MatrixElementCache::setDefectDB)
    .method("setMethod", &MatrixElementCache::setMethod)
    .method("loadElectricDipoleDB", &MatrixElementCache::loadElectricDipoleDB)
    .method("size", &MatrixElementCache::size);

  pi.add_type<StateOne>("StateOne")
    .constructor<std::string, int, int, float, float>()
    .constructor<std::string>()
    .method("string", static_cast<std::string (StateOne::*)() const>(&StateOne::str))
    .method("getN", &StateOne::getN)
    .method("getL", &StateOne::getL)
    .method("getJ", &StateOne::getJ)
    .method("getM", &StateOne::getM)
    .method("getS", &StateOne::getS)
    .method("getSpecies", &StateOne::getSpecies)
    .method("getElement", &StateOne::getElement)
    .method("getEnergy", &StateOne::getEnergy)
    .method("getNStar", &StateOne::getNStar)
    .method("getLabel", &StateOne::getLabel)
    .method("isArtificial", &StateOne::isArtificial)
    .method("isGeneralized", &StateOne::isGeneralized)
    .method("getHash", &StateOne::getHash)
    .method("getReflected", &StateOne::getReflected)
    .method("==", &StateOne::operator==)
    .method("^", &StateOne::operator^)
    .method("!=", &StateOne::operator!=)
    .method("<", &StateOne::operator<)
    .method("<=", &StateOne::operator<=);

  pi.add_type<StateTwo>("StateTwo")
    .constructor<StateOne, StateOne>()
    .method("StateTwo", [](jlcxx::ArrayRef<std::string> init_arr){
      std::array<std::string, 2> str_arr = {init_arr[0], init_arr[1]};
      StateTwo s = StateTwo(str_arr);
      return s;
      })
    .method("StateTwo", [](jlcxx::ArrayRef<std::string> species, 
                           jlcxx::ArrayRef<int> ns, jlcxx::ArrayRef<int> ls, 
                           jlcxx::ArrayRef<float> js, jlcxx::ArrayRef<float> ms){
                             std::array<std::string, 2> species_arr = {species[0], species[1]};
                             std::array<int, 2> n_arr = {ns[0], ns[1]};
                             std::array<int, 2> l_arr = {ls[0], ls[1]};
                             std::array<float, 2> j_arr = {js[0], js[1]};
                             std::array<float, 2> m_arr = {ms[0], ms[1]};
                             StateTwo s = StateTwo(species_arr, n_arr, l_arr, j_arr, m_arr);
                             return s;
                           })
    .method("str", &StateTwo::str)
    .method("getN", [](StateTwo& s){
      std::array<int, 2> n_arr = s.getN();
      return std::make_tuple(n_arr[0], n_arr[1]);
    })
    .method("getN", static_cast<const int& (StateTwo::*)(int) const>(&StateTwo::getN))
    .method("getL", [](StateTwo& s){
      std::array<int, 2> l_arr = s.getL();
      return std::make_tuple(l_arr[0], l_arr[1]);
    })
    .method("getL", static_cast<const int& (StateTwo::*)(int) const>(&StateTwo::getL))
    .method("getJ", [](StateTwo& s){
      std::array<float, 2> j_arr = s.getJ();
      return std::make_tuple(j_arr[0], j_arr[1]);
    })
    .method("getJ", static_cast<const float& (StateTwo::*)(int) const>(&StateTwo::getJ))
    .method("getM", [](StateTwo& s){
      std::array<float, 2> m_arr = s.getM();
      return std::make_tuple(m_arr[0], m_arr[1]);
    })
    .method("getM", static_cast<const float& (StateTwo::*)(int) const>(&StateTwo::getM))
    .method("getS", [](StateTwo& s){
      std::array<float, 2> s_arr = s.getS();
      return std::make_tuple(s_arr[0], s_arr[1]);
    })
    .method("getS", static_cast<const float& (StateTwo::*)(int) const>(&StateTwo::getS))
    .method("getSpecies", [](StateTwo& s){
      std::array<std::string, 2> species_arr = s.getSpecies();
      return std::make_tuple(species_arr[0], species_arr[1]);
    })
    .method("getSpecies", static_cast<const std::string& (StateTwo::*)(int) const>(&StateTwo::getSpecies))
    .method("getElement", [](StateTwo& s){
      std::array<std::string, 2> el_arr = s.getElement();
      return std::make_tuple(el_arr[0], el_arr[1]);
    })
    .method("getElement", static_cast<const std::string& (StateTwo::*)(int) const>(&StateTwo::getElement)) 
    .method("getEnergy", static_cast<double (StateTwo::*)() const>(&StateTwo::getEnergy))
    .method("getEnergy", static_cast<double (StateTwo::*)(int) const>(&StateTwo::getEnergy))
    .method("getNStar", [](StateTwo& s){
      std::array<double, 2> ns_arr = s.getNStar();
      return std::make_tuple(ns_arr[0], ns_arr[1]);
    })
    .method("getNStar", static_cast<double (StateTwo::*)(int) const>(&StateTwo::getNStar))
    .method("getLabel", [](StateTwo& s){
      std::array<std::string, 2> l_arr = s.getLabel();
      return std::make_tuple(l_arr[0], l_arr[1]);
    })
    .method("getLabel", static_cast<const std::string& (StateTwo::*)(int) const>(&StateTwo::getLabel))
    .method("isArtificial", [](StateTwo& s){
      std::array<bool, 2> b_arr = s.isArtificial();
      return std::make_tuple(b_arr[0], b_arr[1]);
    })
    .method("isArtificial", static_cast<bool (StateTwo::*)(int) const>(&StateTwo::isArtificial))
    .method("isGeneralized", [](StateTwo& s){
      std::array<bool, 2> b_arr = s.isGeneralized();
      return std::make_tuple(b_arr[0], b_arr[1]);
    })
    .method("isGeneralized", static_cast<bool (StateTwo::*)(int) const>(&StateTwo::isGeneralized))
    .method("getFirstState", &StateTwo::getFirstState)
    .method("getSecondState", &StateTwo::getSecondState)
    .method("getHash", &StateTwo::getHash)
    .method("getReflected", &StateTwo::getReflected)
    .method("==", &StateTwo::operator==)
    .method("^", &StateTwo::operator^)
    .method("!=", &StateTwo::operator!=)
    .method("<", &StateTwo::operator<)
    .method("<=", &StateTwo::operator<=);

  pi.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("SystemBase")
    .apply<SystemBase<StateOne>>([](auto wrapped){
      typedef typename decltype(wrapped)::type WrappedT;
    })
    .apply<SystemBase<StateTwo>>([](auto wrapped){
      typedef typename decltype(wrapped)::type WrappedT;
    });

  pi.add_type<eigen_sparse_t>("eigen_sparse_t")
    .method("nonzerorealvalues", [](eigen_sparse_t& e){
      jlcxx::Array<double> ret;
      #ifdef USE_COMPLEX
      for (int i=0; i<e.nonZeros(); i++){
          ret.push_back((e.valuePtr()[i]).real());
      }
      #else
      for (int i=0; i<e.nonZeros(); i++){
          ret.push_back(e.valuePtr()[i]);
      }
      #endif     
      return ret;
    })
    .method("nonzeroimagvalues", [](eigen_sparse_t& e){
      jlcxx::Array<double> ret;
      #ifdef USE_COMPLEX
      for (int i=0; i<e.nonZeros(); i++){
          ret.push_back((e.valuePtr()[i]).imag());
      }
      #else
      for (int i=0; i<e.nonZeros(); i++){
          ret.push_back(0);
      }
      #endif     
      return ret;
    })
    .method("outerIndex", [](eigen_sparse_t& e){
      jlcxx::Array<int> ret;
      for (int i=0; i<e.outerSize(); i++){
        ret.push_back(e.outerIndexPtr()[i]);
      }
      return ret;
    })
    .method("innerIndex", [](eigen_sparse_t& e){
      jlcxx::Array<int> ret;
      for (int i=0; i<e.nonZeros(); i++){
        ret.push_back(e.innerIndexPtr()[i]);
      }
      return ret;
    });

  pi.add_type<SystemOne>("SystemOne", jlcxx::julia_type<SystemBase<StateOne>>())
    .constructor<std::string, MatrixElementCache&>()
    .constructor<std::string, MatrixElementCache&, bool>()
    // SystemBase methods ///////////////////////////////////////////////////////////////////////
    .method("restrictEnergy", &SystemOne::restrictEnergy)
    .method("restrictN", static_cast<void (SystemOne::*)(int, int)>(&SystemOne::restrictN))
    .method("restrictN", [](SystemOne& s, jlcxx::ArrayRef<int> n_jl){
      std::set<int> n = {n_jl[0], n_jl[1]};
      s.restrictN(n);
    })
    .method("restrictL", static_cast<void (SystemOne::*)(int, int)>(&SystemOne::restrictL))
    .method("restrictL", [](SystemOne& s, jlcxx::ArrayRef<int> l_jl){
      std::set<int> l = {l_jl[0], l_jl[1]};
      s.restrictL(l);
    })
    .method("restrictJ", static_cast<void (SystemOne::*)(float, float)>(&SystemOne::restrictJ))
    .method("restrictJ", [](SystemOne& s, jlcxx::ArrayRef<float> j_jl){
      std::set<float> j = {j_jl[0], j_jl[1]};
      s.restrictJ(j);
    })
    .method("restrictM", static_cast<void (SystemOne::*)(float, float)>(&SystemOne::restrictM))
    .method("restrictM", [](SystemOne& s, jlcxx::ArrayRef<float> m_jl){
      std::set<float> m = {m_jl[0], m_jl[1]};
      s.restrictM(m);
    })
    .method("diagonalize", static_cast<void (SystemOne::*)()>(&SystemOne::diagonalize))
    .method("getHamiltonian", static_cast<eigen_sparse_t& (SystemOne::*)()>(&SystemOne::getHamiltonian))
    /////////////////////////////////////////////////////////////////////////////////////////////
    .method("getSpecies", &SystemOne::getSpecies)
    .method("setEfield", [](SystemOne& s, jlcxx::ArrayRef<double> efield){
      std::array<double, 3> field = {efield[0], efield[1], efield[2]};
      s.setEfield(field);
    })
    .method("setEfield", [](SystemOne& s, jlcxx::ArrayRef<double> efield, 
                                          jlcxx::ArrayRef<double> z_axis, 
                                          jlcxx::ArrayRef<double> y_axis){
      std::array<double, 3> field = {efield[0], efield[1], efield[2]};
      std::array<double, 3> to_z_axis = {z_axis[0], z_axis[1], z_axis[2]};
      std::array<double, 3> to_y_axis = {y_axis[0], y_axis[1], y_axis[2]};
      s.setEfield(field, to_z_axis, to_y_axis);
    })
    .method("setEfield", [](SystemOne& s, jlcxx::ArrayRef<double> efield, double alpha, double beta, double gamma){
      std::array<double, 3> field = {efield[0], efield[1], efield[2]};
      s.setEfield(field, alpha, beta, gamma);
    })
    .method("setBfield", [](SystemOne& s, jlcxx::ArrayRef<double> bfield){
      std::array<double, 3> field = {bfield[0], bfield[1], bfield[2]};
      s.setBfield(field);
    })
    .method("setBfield", [](SystemOne& s, jlcxx::ArrayRef<double> bfield, 
                                          jlcxx::ArrayRef<double> z_axis, 
                                          jlcxx::ArrayRef<double> y_axis){
      std::array<double, 3> field = {bfield[0], bfield[1], bfield[2]};
      std::array<double, 3> to_z_axis = {z_axis[0], z_axis[1], z_axis[2]};
      std::array<double, 3> to_y_axis = {y_axis[0], y_axis[1], y_axis[2]};
      s.setBfield(field, to_z_axis, to_y_axis);
    })
    .method("setBfield", [](SystemOne& s, jlcxx::ArrayRef<double> bfield, double alpha, double beta, double gamma){
      std::array<double, 3> field = {bfield[0], bfield[1], bfield[2]};
      s.setBfield(field, alpha, beta, gamma);
    })
    .method("enableDiamagnetism", &SystemOne::enableDiamagnetism)
    .method("setConservedParityUnderReflection", &SystemOne::setConservedParityUnderReflection)
    .method("setConservedMomentaUnderRotation", [](SystemOne& s, jlcxx::ArrayRef<float> momenta_jl){
      std::set<float> momenta;
      for (unsigned i=0; i<momenta_jl.size(); i++){
        momenta.insert(momenta_jl[i]);
      }
      const std::set<float> &momenta_const = momenta;
      s.setConservedMomentaUnderRotation(momenta_const);
    });

  pi.add_type<SystemTwo>("SystemTwo", jlcxx::julia_type<SystemBase<StateTwo>>())
    .constructor<SystemOne, SystemOne, MatrixElementCache&>()
    .constructor<SystemOne, SystemOne, MatrixElementCache&, bool>()
    .constructor<SystemTwo>()
    // SystemBase methods ///////////////////////////////////////////////////////////////////////
    .method("restrictEnergy", &SystemTwo::restrictEnergy)
    .method("restrictN", static_cast<void (SystemTwo::*)(int, int)>(&SystemTwo::restrictN))
    .method("restrictN", [](SystemTwo& s, jlcxx::ArrayRef<int> n_jl){
      std::set<int> n = {n_jl[0], n_jl[1]};
      s.restrictN(n);
    })
    .method("restrictL", static_cast<void (SystemTwo::*)(int, int)>(&SystemTwo::restrictL))
    .method("restrictL", [](SystemTwo& s, jlcxx::ArrayRef<int> l_jl){
      std::set<int> l = {l_jl[0], l_jl[1]};
      s.restrictL(l);
    })
    .method("restrictJ", static_cast<void (SystemTwo::*)(float, float)>(&SystemTwo::restrictJ))
    .method("restrictJ", [](SystemTwo& s, jlcxx::ArrayRef<float> j_jl){
      std::set<float> j = {j_jl[0], j_jl[1]};
      s.restrictJ(j);
    })
    .method("restrictM", static_cast<void (SystemTwo::*)(float, float)>(&SystemTwo::restrictM))
    .method("restrictM", [](SystemTwo& s, jlcxx::ArrayRef<float> m_jl){
      std::set<float> m = {m_jl[0], m_jl[1]};
      s.restrictM(m);
    })
    .method("diagonalize", static_cast<void (SystemTwo::*)()>(&SystemTwo::diagonalize))
    .method("getHamiltonian", static_cast<eigen_sparse_t& (SystemTwo::*)()>(&SystemTwo::getHamiltonian))
    .method("getOverlap", [](SystemTwo& st, StateTwo& s, double alpha, double beta, double gamma){
      eigen_vector_double_t overlap = st.getOverlap(s, alpha, beta, gamma);
      jlcxx::Array<double> ret;
      for (unsigned i=0; i<overlap.size(); i++){
        ret.push_back(overlap.data()[i]);
      }
      return ret;
    })
    /////////////////////////////////////////////////////////////////////////////////////////////
    .method("getSpecies", [](SystemTwo& s){
      std::array<std::string, 2> species_arr = s.getSpecies();
      return std::make_tuple(species_arr[0], species_arr[1]);
    })
    .method("getStatesFirst", &SystemTwo::getStatesFirst)
    .method("getStatesSecond", &SystemTwo::getStatesSecond)
    .method("enableGreenTensor", &SystemTwo::enableGreenTensor)
    .method("setSurfaceDistance", &SystemTwo::setSurfaceDistance)
    .method("setAngle", &SystemTwo::setAngle)
    .method("setDistance", &SystemTwo::setDistance)
    .method("setDistanceVector", [](SystemTwo& s, jlcxx::ArrayRef<double> dvec){
      std::array<double, 3> d = {dvec[0], dvec[1], dvec[2]};
      s.setDistanceVector(d);
    })
    .method("setOrder", &SystemTwo::setOrder)
    .method("setConservedParityUnderPermutation", &SystemTwo::setConservedParityUnderPermutation)
    .method("setConservedParityUnderInversion", &SystemTwo::setConservedParityUnderInversion)
    .method("setConservedParityUnderReflection", &SystemTwo::setConservedParityUnderReflection)
    .method("setConservedMomentaUnderRotation", [](SystemTwo& s, jlcxx::ArrayRef<int> momenta_jl){
      std::set<int> momenta;
      for (unsigned i=0; i<momenta_jl.size(); i++){
        momenta.insert(momenta_jl[i]);
      }
      const std::set<int> &momenta_const = momenta;
      s.setConservedMomentaUnderRotation(momenta_const);
    });

  pi.add_type<QuantumDefect>("QuantumDefect")
    .constructor<std::string const&, int, int, double>()
    .constructor<std::string const&, int, int, double, std::string const&>()
    .method("n", [](QuantumDefect &qd){return qd.n;})
    .method("l", [](QuantumDefect &qd){return qd.l;})
    .method("j", [](QuantumDefect &qd){return qd.j;})
    .method("ac", [](QuantumDefect &qd){return qd.ac;})
    .method("Z", [](QuantumDefect &qd){return qd.Z;})
    .method("a1", [](QuantumDefect &qd){return qd.a1;})
    .method("a2", [](QuantumDefect &qd){return qd.a2;})
    .method("a3", [](QuantumDefect &qd){return qd.a3;})
    .method("a4", [](QuantumDefect &qd){return qd.a4;})
    .method("rc", [](QuantumDefect &qd){return qd.rc;})
    .method("nstar", [](QuantumDefect &qd){return qd.nstar;})
    .method("energy", [](QuantumDefect &qd){return qd.energy;});
}
