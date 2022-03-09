/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * NonlinearStrainOperatorFast.h
 *
 * Code generation for function 'NonlinearStrainOperatorFast'
 *
 */

#ifndef NONLINEARSTRAINOPERATORFAST_H
#define NONLINEARSTRAINOPERATORFAST_H

/* Include files */
#include "NonlinearStrainOperatorFast_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void NonlinearStrainOperatorFast(const emxArray_real_T *N,
                                        const emxArray_real_T *dNdx,
                                        const double F[9], emxArray_real_T *Bn,
                                        emxArray_real_T *Bg,
                                        emxArray_real_T *NN, double *tau);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (NonlinearStrainOperatorFast.h) */
