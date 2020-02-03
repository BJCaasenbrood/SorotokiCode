/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_minv_api.h
 *
 * Code generation for function '_coder_minv_api'
 *
 */

#ifndef _CODER_MINV_API_H
#define _CODER_MINV_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_minv_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void minv(real_T A[9], real_T X[9]);
extern void minv_api(const mxArray * const prhs[1], int32_T nlhs, const mxArray *
                     plhs[1]);
extern void minv_atexit(void);
extern void minv_initialize(void);
extern void minv_terminate(void);
extern void minv_xil_terminate(void);

#endif

/* End of code generation (_coder_minv_api.h) */
