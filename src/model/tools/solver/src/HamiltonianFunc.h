#ifndef HAMILTONIAN_FUNC
#define HAMILTONIAN_FUNC

#include <Eigen/Dense>
#include "AutoDiffScalar.h"
#include "AutoDiffJacobian.h"

template <short N_IN>
class Hamiltonian {
   
  public:
   
    enum
    {
      InputsAtCompileTime = N_IN,
      ValuesAtCompileTime = 1
    };

    typedef Eigen::Matrix<double, N_IN, 1> InputType;
    typedef Eigen::Matrix<double, 1, 1> ValueType;
   
    template <typename T1, typename T2>
    void operator()
    (
      const Eigen::Matrix<T1, N_IN, 1>& x, 
      Eigen::Matrix<T2, 1, 1>* H
    ) const
    {

      Eigen::Matrix<T1, N_IN/2, 1> q;
      Eigen::Matrix<T1, N_IN/2, 1> p;

      q.noalias() = x.block(0,0,N_IN/2,1);
      p.noalias() = x.block(N_IN/2,0,N_IN/2,1);

      (*H) << 0.5*q.transpose()*q + 0.5*p.transpose()*p;
    } 

};

#endif