/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * computeForwardKinematicsFast.h
 *
 * Code generation for function 'computeForwardKinematicsFast'
 *
 */

#ifndef COMPUTEFORWARDKINEMATICSFAST_H
#define COMPUTEFORWARDKINEMATICSFAST_H

/* Include files */
#include "computeForwardKinematicsFast_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void computeForwardKinematicsFast(
    const emxArray_real_T *x, double ds, const double p0[3],
    const double Phi0[9], const emxArray_real_T *xia0,
    const emxArray_real_T *Th, const emxArray_real_T *Ba, emxArray_real_T *gtmp,
    emxArray_real_T *Jtmp);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (computeForwardKinematicsFast.h) */
