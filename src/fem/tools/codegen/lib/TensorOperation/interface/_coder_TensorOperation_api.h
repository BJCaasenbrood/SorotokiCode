/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_TensorOperation_api.h
 *
 * Code generation for function '_coder_TensorOperation_api'
 *
 */

#ifndef _CODER_TENSOROPERATION_API_H
#define _CODER_TENSOROPERATION_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_TensorOperation_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void TensorOperation(real_T A[9], real_T B[9], boolean_T Arg, real_T T[36]);
extern void TensorOperation_api(const mxArray * const prhs[3], int32_T nlhs,
  const mxArray *plhs[1]);
extern void TensorOperation_atexit(void);
extern void TensorOperation_initialize(void);
extern void TensorOperation_terminate(void);
extern void TensorOperation_xil_terminate(void);

#endif

/* End of code generation (_coder_TensorOperation_api.h) */
