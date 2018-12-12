
#ifndef greentensor_h
#define greentensor_h

// #include "SystemTwo.h"
// #include "dtypes.h"

#include<math.h>
#include<stdio.h>
#include<fstream>//um Datei zu erzeugen, zu Testzwecken
#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<ios>
#include<cstdlib>
#include<string>
#include<cmath>
#include<complex>
#include <Eigen/Sparse>
#include <vector>


class GreenTensor{

    public:
    Eigen::Matrix<std::complex<double>, 3, 3> tensor;
    
    void vacuum(double x,double z);
    void plate(double x, double zA, double zB );
    
    GreenTensor(double x, double z);     
    GreenTensor(double x, double zA, double zB);
};

void GreenTensor::vacuum(double x, double z){
      double distance = std::sqrt(x*x + z*z);
      double Kdelta;
      double vecrho[3];
      z = std::abs(z);
      vecrho[0] = x/(distance*distance);
      vecrho[1] = 0.;
      vecrho[2] = z/(distance*distance);
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  Kdelta = 0.;
	  if(i==j){
	    Kdelta = 1.;
	  }
	  tensor(i,j) = (Kdelta - 3.*vecrho[i]*vecrho[j])/std::pow(distance,3.);
	}
      }
}

void GreenTensor::plate(double x, double zA, double zB){
      if(zA<0. || zB<0.){
	std::cout<<"error! z<0"<<std::endl;
      }
      double zp = zA+zB;
      double rp = std::sqrt(x*x+zp*zp);
      tensor(0,0) += (1. - (3.*x*x)/(rp*rp))*std::pow(rp,-3.);
      tensor(0,2) +=  ((3.*x*zp)/(rp*rp))*std::pow(rp,-3.);
      tensor(1,1) += 1.*std::pow(rp,-3.);
      tensor(2,0) += (- 3.*x*zp/(rp*rp))*std::pow(rp,-3.);
      tensor(2,2) += (2. - 3.*x*x/(rp*rp))*std::pow(rp,-3.);
}
GreenTensor::GreenTensor(double x, double z){
      for(int i =0;i<3;i++){
	for(int j = 0;j<3;j++){
	  tensor(i,j) = 0.;
	}
      }
      vacuum(x,z);
}

GreenTensor::GreenTensor(double x, double zA, double zB){
      for(int i =0;i<3;i++){
	for(int j = 0;j<3;j++){
	  tensor(i,j) = 0.;
	}
      }
      vacuum(x,zA-zB);
      plate(x,zA,zB);      
}

#endif
