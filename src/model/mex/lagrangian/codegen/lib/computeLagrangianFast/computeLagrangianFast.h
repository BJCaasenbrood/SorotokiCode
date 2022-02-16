/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * computeLagrangianFast.h
 *
 * Code generation for function 'computeLagrangianFast'
 *
 */

#ifndef COMPUTELAGRANGIANFAST_H
#define COMPUTELAGRANGIANFAST_H

/* Include files */
#include "computeLagrangianFast_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void computeLagrangianFast(
    const emxArray_real_T *x, const emxArray_real_T *dx, double ds,
    const double p0[3], const double Phi0[9], const emxArray_real_T *xia0,
    const emxArray_real_T *Th, const emxArray_real_T *Ba, const double Ktt[36],
    const double Mtt[36], double Zeta, emxArray_real_T *M, emxArray_real_T *C,
    emxArray_real_T *K, emxArray_real_T *R, emxArray_real_T *G, double p[3],
    double Phi[9], emxArray_real_T *J, double *Vg, double *Kin);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (computeLagrangianFast.h) */
