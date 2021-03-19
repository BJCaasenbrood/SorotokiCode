/*
// solvers for Algebraic Riccati equation
// - Iteration (continuous)
// - Iteration (discrete)
// - Arimoto-Potter
//
// author: Horibe Takamasa
*/

#include "CARE.h"

bool solveRiccatiIterationC(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                            const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                            Eigen::MatrixXd &P, const double dt,
                            const double &tolerance,
                            const uint iter_max) {
  P = Q; // initialize

  Eigen::MatrixXd P_next;

  Eigen::MatrixXd AT = A.transpose();
  Eigen::MatrixXd BT = B.transpose();
  Eigen::MatrixXd Rinv = R.inverse();

  double diff;
  for (uint i = 0; i < iter_max; ++i) {
    P_next = P + (P * A + AT * P - P * B * Rinv * BT * P + Q) * dt;
    diff = fabs((P_next - P).maxCoeff());
    P = P_next;
    if (diff < tolerance) {
      //std::cout << "iteration mumber = " << i << std::endl;
      return true;
    }

  }
  return false; // over iteration limit
}

bool solveRiccatiIterationD(const Eigen::MatrixXd &Ad,
                            const Eigen::MatrixXd &Bd, const Eigen::MatrixXd &Q,
                            const Eigen::MatrixXd &R, Eigen::MatrixXd &P,
                            const double &tolerance,
                            const uint iter_max) {
  P = Q; // initialize

  Eigen::MatrixXd P_next;

  Eigen::MatrixXd AdT = Ad.transpose();
  Eigen::MatrixXd BdT = Bd.transpose();
  Eigen::MatrixXd Rinv = R.inverse();

  double diff;
  for (uint i = 0; i < iter_max; ++i) {
    // -- discrete solver --
    P_next = AdT * P * Ad -
             AdT * P * Bd * (R + BdT * P * Bd).inverse() * BdT * P * Ad + Q;

    diff = fabs((P_next - P).maxCoeff());
    P = P_next;
    if (diff < tolerance) {
      //std::cout << "iteration mumber = " << i << std::endl;
      return true;
    }
  }
  return false; // over iteration limit
}

void solveRiccatiArimotoPotter(Eigen::MatrixXd A,
                               Eigen::MatrixXd B,
                               Eigen::MatrixXd Q,
                               Eigen::MatrixXd R, Eigen::MatrixXd &P) {

  uint dim_x = A.rows();
  uint dim_u = B.cols();

  // set Hamilton matrix
  Eigen::MatrixXd Ham(2 * dim_x, 2 * dim_x);
  Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
  Ham << A, -(B * R.inverse())* B.transpose(), -Q, -A.transpose();

  // calc eigenvalues and eigenvectors
  Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

  // check eigen values
  // std::cout << "eigen values：\n" << Eigs.eigenvalues() << std::endl;
  // std::cout << "eigen vectors：\n" << Eigs.eigenvectors() << std::endl;

  // extract stable eigenvectors into 'eigvec'
  Eigen::MatrixXcd eigvec(2 * dim_x, dim_x);
  eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);

  int j = 0;
  for (int i = 0; i < 2 * dim_x; ++i) {
    if (Eigs.eigenvalues()[i].real() < -1.E-5) { // was 0.
      eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
      ++j;
    }
  }

  // calc P with stable eigen vector matrix
  Eigen::MatrixXcd Vs_1(dim_x,dim_x), Vs_2(dim_x,dim_x);
  Vs_1.noalias() = eigvec.block(0, 0, dim_x, dim_x);
  Vs_2.noalias()  = eigvec.block(dim_x, 0, dim_x, dim_x);
  P.noalias() = (Vs_2 * Vs_1.inverse()).real();

  //return true;
}
