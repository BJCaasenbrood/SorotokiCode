/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PolarDecompositionFast_api.h
 *
 * Code generation for function 'PolarDecompositionFast'
 *
 */

#ifndef _CODER_POLARDECOMPOSITIONFAST_API_H
#define _CODER_POLARDECOMPOSITIONFAST_API_H

/* Include files */
#include "emlrt.h"
#include "tmwtypes.h"
#include <string.h>

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void PolarDecompositionFast(real_T F[9], creal_T R[9], creal_T S[9],
                            creal_T V[9]);

void PolarDecompositionFast_api(const mxArray *prhs, int32_T nlhs,
                                const mxArray *plhs[3]);

void PolarDecompositionFast_atexit(void);

void PolarDecompositionFast_initialize(void);

void PolarDecompositionFast_terminate(void);

void PolarDecompositionFast_xil_shutdown(void);

void PolarDecompositionFast_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_PolarDecompositionFast_api.h) */
