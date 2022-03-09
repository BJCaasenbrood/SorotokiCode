/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PiollaStressNHFast.h
 *
 * Code generation for function 'PiollaStressNHFast'
 *
 */

#ifndef PIOLLASTRESSNHFAST_H
#define PIOLLASTRESSNHFAST_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void PiollaStressNHFast(const double F[9], double C100, double K,
                               double S[9], double D[36], double *P);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (PiollaStressNHFast.h) */
