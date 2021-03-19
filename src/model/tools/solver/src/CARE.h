/*
// solvers for Algebraic Riccati equation
// - Iteration (continuous)
// - Iteration (discrete)
// - Arimoto-Potter
//
// author: Horibe Takamasa
*/
//#pragma once
#ifndef RICCATI_H
#define RICCATI_H


#include <Eigen/Dense>
#include <iostream>
#include <time.h>
#include <vector>

/* Itereation method for continuous model */
bool solveRiccatiIterationC(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                            const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                            Eigen::MatrixXd &P, const double dt = 0.001,
                            const double &tolerance = 1.E-5,
                            const uint iter_max = 100000);

/* Itereation method for discrete model */
bool solveRiccatiIterationD(const Eigen::MatrixXd &Ad,
                            const Eigen::MatrixXd &Bd, const Eigen::MatrixXd &Q,
                            const Eigen::MatrixXd &R, Eigen::MatrixXd &P,
                            const double &tolerance = 1.E-5,
                            const uint iter_max = 100000);

/* Arimoto Potter method for continuous model */
void solveRiccatiArimotoPotter(Eigen::MatrixXd A,
                               Eigen::MatrixXd B,
                               Eigen::MatrixXd Q,
                               Eigen::MatrixXd R, Eigen::MatrixXd &P);

#endif