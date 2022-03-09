/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PolarDecompositionFast_mex.h
 *
 * Code generation for function 'PolarDecompositionFast'
 *
 */

#ifndef _CODER_POLARDECOMPOSITIONFAST_MEX_H
#define _CODER_POLARDECOMPOSITIONFAST_MEX_H

/* Include files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS(void);

void unsafe_PolarDecompositionFast_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                               int32_T nrhs,
                                               const mxArray *prhs[1]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_PolarDecompositionFast_mex.h) */
