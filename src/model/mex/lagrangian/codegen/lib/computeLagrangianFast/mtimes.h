/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.h
 *
 * Code generation for function 'mtimes'
 *
 */

#ifndef MTIMES_H
#define MTIMES_H

/* Include files */
#include "computeLagrangianFast_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_mtimes(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *C);

void mtimes(const emxArray_real_T *A, const double B[36], emxArray_real_T *C);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (mtimes.h) */
