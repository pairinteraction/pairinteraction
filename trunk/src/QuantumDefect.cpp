#include "QuantumDefect.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <exception>


struct no_defect : public std::exception {
  const char* what () const throw ()
  {
    return "There is no defect available";
  }
};


void QuantumDefect::Rb87(int n, int l, double j) {
  //***************
  //* Rubidium 87 *
  //***************
  ac_ = 9.076;
  Z_ = 37;
  switch (l) {
  case 0: // S-state
    a1_ = 3.69628474;
    a2_ = 1.64915255;
    a3_ = -9.86069196;
    a4_ = 0.19579987;
    rc_ = 1.66242117;
    break;
  case 1: // P-state
    a1_ = 4.44088978;
    a2_ = 1.92828831;
    a3_ = -16.79597770;
    a4_ = -0.81633314;
    rc_ = 1.50195124;
    break;
  case 2: // D-state
    a1_ = 3.78717363;
    a2_ = 1.57027864;
    a3_ = -11.6558897;
    a4_ = 0.52942835;
    rc_ = 4.86851938;
    break;
  default: // L>= 3
    a1_ = 2.39848933;
    a2_ = 1.76810544;
    a3_ = -12.0710678;
    a4_ = 0.77256589;
    rc_ = 4.79831327;
    break;
  }

  if ( n == 5 ) {
    switch (l) {
    case 0:
      if ( j == 0.5 ) d0 = 3.1944;
      break;
    case 1:
      if ( j == 0.5 ) d0 = 2.71906;
      else if ( j == 1.5 ) d0 = 2.7071;
      break;
    }
  }
  else if ( n == 6 ) {
    switch (l) {
    case 1:
      if ( j == 1.5 ) d0 = 2.66885;
      break;
    }
  }
  else if ( n > 19 ) {
    switch (l) {
    case 0:
      d0 = 3.1311807;
      d2 = 0.1787;
      break;
    case 1:
      if ( j == 0.5 ) {
        d0 = 2.6548849;
        d2 = 0.29;
      }
      else if ( j == 1.5 ) {
        d0 = 2.6416737;
        d2 = 0.295;
      }
      break;
    case 2:
      if ( j == 1.5 ) {
        d0 =  1.34809171;
        d2 = -0.60286;
      }
      if ( j == 2.5 ) {
        d0 = 1.34646572;
        d2 = -0.596;
      }
      break;
    case 3:
      if ( j == 2.5 ) {
        d0 = 0.0165192;
        d2 = -0.085;
      }
      if ( j == 3.5 ) {
        d0 = 0.0165437;
        d2 = -0.086;
      }
      break;
    default:
      d0 = 0;
      break;
    }
  }    
}


void QuantumDefect::Rb85(int n, int l, double j) {
  //***************
  //* Rubidium 85 *
  //***************
  ac_ = 9.076;
  Z_ = 35;        
  switch (l) {
  case 0: // S-state
    a1_ = 3.69628474;
    a2_ = 1.64915255;
    a3_ = -9.86069196;
    a4_ = 0.19579987;
    rc_ = 1.66242117;
    break;
  case 1: // P-state
    a1_ = 4.44088978;
    a2_ = 1.92828831;
    a3_ = -16.79597770;
    a4_ = -0.81633314;
    rc_ = 1.50195124;
    break;
  case 2: // D-state
    a1_ = 3.78717363;
    a2_ = 1.57027864;
    a3_ = -11.6558897;
    a4_ = 0.52942835;
    rc_ = 4.86851938;
    break;
  default: //  L>= 3
    a1_ = 2.39848933;
    a2_ = 1.76810544;
    a3_ = -12.0710678;
    a4_ = 0.77256589;
    rc_ = 4.79831327;
    break;
  }
}


void QuantumDefect::Cs(int n, int l, double j) {
  //***********
  //* Caesium *
  //***********
  ac_ = 15.644;
  Z_ = 55;        
  switch (l) {
  case 0: // S-state
    a1_ = 3.49546309;
    a2_ = 1.475338;
    a3_ = -9.72143084;
    a4_ = 0.02629242;
    rc_ = 1.9204693;
    break;
  case 1: // P-state
    a1_ = 4.69366096;
    a2_ = 1.71398344;
    a3_ = -24.6562428;
    a4_ = -0.09543125;
    rc_ = 2.13383095;
    break;
  case 2: // D-state
    a1_ = 4.32466196;
    a2_ = 1.61365288;
    a3_ = -6.7012885;
    a4_ = -0.74095193;
    rc_ = 0.93007296;
    break;
  default: //  L>= 3
    a1_ = 3.01048361;
    a2_ = 1.40000001;
    a3_ = -3.20036138;
    a4_ = 0.00034538;
    rc_ = 1.99969677;
    break;
  }

  if (n==5) {
    // from  J. Phys. Chem. Ref. Data 38, 761 (2009); doi: 10.1063/1.3132702 
    // with Rydberg constant and ionization energy from Niemax and Lorenzen
    switch (l) {
    case 2:
      if ( j == 1.5 ) { // Cs133 5D_3/2 state
        d0 = 2.4523471;
        d2 = 0;
      }
      else if ( j == 2.5 ) { //Cs133 5D_5/2 state
        d0 = 2.44496294;
        d2 = 0;
      }
      break;
    }
  }
  else if ( n >= 6 ) {
    switch (l) {
    case 0: // nS Rydberg states (from Niemax Lorenzen, doi: 10.1007/BF01419370)
      if ( j == 0.5 ) {
        d0 = 4.0493527;
        d2 = 0.2381;
        d4 = 0.24688;
        d6 = 0.06785;
        d8 = 0.1135; 
      }
      break;
    case 1: // nP Rydberg states (from Niemax Lorenzen, doi: 10.1007/BF01419370)
      if ( j == 0.5 ) { //nP_1/2 states
        d0 = 3.5914856;
        d2 = 0.380223; 
        d4 = -0.64771; 
        d6 = 20.9538; 
        d8 = -84.108; 
      }
      else if ( j == 1.5 ) {//nP_3/2 states
        d0 = 3.5589599;
        d2 = 0.392469; 
        d4 = -0.67431; 
        d6 = 22.3531; 
        d8 = -92.289; 
      }
      break;
    case 2: //nD Rydberg states (from Niemax Lorenzen, doi: 10.1007/BF01419370)
      if ( j == 1.5 ) { // nD_3/2 states
        d0 = 2.4754562;
        d2 = 0.00932;  
        d4 = -0.43498; 
        d6 = -0.76358; 
        d8 = -18.0061; 
      }
      else if ( j == 2.5 ) {// nD_5/2 states
        d0 = 2.4663091;
        d2 = 0.014964; 
        d4 = -0.45828; 
        d6 = -0.25489; 
        d8 = -19.6900; 
      }
      break;
    case 3: //nF Rydberg states (from Goy and Haroche PRA 26, 27332742 (1982))
      if ( j == 2.5 ) { // nF_5/2 states
        d0 = 0.033392;
        d2 = -0.00191;
      }
      else if ( j == 3.5 ) { // nF_7/2 states
        d0 = -0.00191;
        d2 = -0.00191; 
      }
      break;
    default: //for L>3
      d0=0;
      d2=0;
      break;
    }
  }

  if ( (n==7) && (l==1) ) {
    // from  J. Phys. Chem. Ref. Data 38, 761 (2009); doi: 10.1063/1.3132702 
    // with Rydberg constant and ionization energy from Niemax and Lorenzen
    if ( j == 0.5 ) { //Cs133 7P_1/2 state
      d0 = 3.626253245;
      d2 = 0;
    }
    else if ( j == 1.5 ) { //Cs133 7P_3/2 state
      d0 = 3.594122558;
      d2 = 0;
    }
  }
  else if ( (n==8) && (l==1) ) {
    // from  J. Phys. Chem. Ref. Data 38, 761 (2009); doi: 10.1063/1.3132702 
    // with Rydberg constant and ionization energy from Niemax and Lorenzen
    if ( j == 0.5 ) { //Cs133 8P_1/2 state
      d0 = 3.611362476;
      d2 = 0;
    }
    if ( j == 1.5 ) { //Cs133 8P_3/2 state
      d0 = 3.57917986;
      d2 = 0;
    }
  }
}


void QuantumDefect::calc_delta(int n) {
  std::stringstream warning;
  bool warn = false;
  double denom;
  warning << "The following coefficients are not given:";
  delta_ = 0;
  if ( d0 > 0 ) {
    delta_ += d0;
    
    denom = (n - d0)*(n-d0);
    if ( d2 >= 0 )
      delta_ += d2 / denom;
    else { warn = true; warning << " d2"; }
    
    denom *= (n - d0)*(n-d0);
    if ( d4 >= 0 )
      delta_ += d4 / denom;
    else { warn = true; warning << " d4"; }
    
    denom *= (n - d0)*(n-d0);
    if ( d6 >= 0 )
      delta_ += d6 / denom;
    else { warn = true; warning << " d6"; }
    
    denom *= (n - d0)*(n-d0);
    if ( d8 >= 0 )
      delta_ += d8 / denom;
    else { warn = true; warning << " d8"; }
  } else throw no_defect();

  if (warn)
    std::cerr << warning.str() << std::endl;
}


QuantumDefect::QuantumDefect(std::string species, int n, int l, double j)
  : d0(-1), d2(-1), d4(-1), d6(-1), d8(-1), delta(delta_),
    ac(ac_), Z(Z_), a1(a1_), a2(a2_), a3(a3_), a4(a4_), rc(rc_)
{
  if (species == std::string("Rb87"))
    Rb87(n, l, j);
  else if (species == std::string("Rb85"))
    Rb85(n, l, j);
  else if (species == std::string("Cs"))
    Cs(n, l, j);
  else throw no_defect();

  calc_delta(n);
}

/*
// resulting potential including finite core size and core polarizibility
Znl _ =  1 + (Z-1) *exp(-a1*x.^2) - x.^2 .* (a3 + a4*x.^2) .*exp(-a2*x.^2);
Vc _ =  - Znl./(x.^2) - (ac./(2*x.^8)) .* (1-exp(-((x.^2)/rc).^6));

// ------------------------------------------------
// for S,P,D,F states finestructure is included
if (L<4)
    alpha_ = 7.29735257e-3;
    Vso_ = alpha.^2./(4* x.^6) *(J*(J+1)-L*(L+1)-0.5*(1.5));
else
    Vso_ = 0;
end

Vx _ =  Vc+Vso; // total potential in the radial Schrödinger equation           
        
end

*/
