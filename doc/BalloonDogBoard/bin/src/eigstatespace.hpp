#ifndef _EIGEN_QUADSOLVE_HPP_
#define _EIGEN_QUADSOLVE_HPP_

#include <Arduino.h>
#include <Dense>
#include <cmath>

namespace Eigen {

inline void statespace(const float t, const VectorXd &Y, const VectorXd &U, VectorXd &Y_p){
	
  const float mu = 0.01;
  
  Y_p(0) = Y(1);
  Y_p(1) = mu*(1-Y(0)*Y(0))*Y(1) - Y(0);
}

inline void update_statespace(const float t, const float dt, VectorXd &Y, const VectorXd &U){

  int N = Y.size();

  VectorXd K1(N);
  VectorXd K2(N);
  VectorXd K3(N);
  VectorXd K4(N);
  VectorXd Ytmp(N);
  
  statespace(t, Y, U, K1);
  Ytmp = Y + 0.5*dt*K1;
  
  statespace(t + 0.5*dt, Ytmp, U, K2);
  Ytmp = Y + 0.5*dt*K2;
  
  statespace(t + 0.5*dt, Ytmp, U, K3);
  Ytmp = Y + 0.5*dt*K3;
  
  statespace(t + dt, Ytmp, U, K4);
  Y = Y + (dt / 6.0) * (K1 + 2.0*(K2 + K3) + K4);
}

}
#endif
