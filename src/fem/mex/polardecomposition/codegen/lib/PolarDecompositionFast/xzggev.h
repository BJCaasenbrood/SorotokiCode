/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzggev.h
 *
 * Code generation for function 'xzggev'
 *
 */

#ifndef XZGGEV_H
#define XZGGEV_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void xzggev(creal_T A[9], int *info, creal_T alpha1[3], creal_T beta1[3],
            creal_T V[9]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xzggev.h) */
