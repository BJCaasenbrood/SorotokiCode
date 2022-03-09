/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PiollaStressNHFast_api.h
 *
 * Code generation for function 'PiollaStressNHFast'
 *
 */

#ifndef _CODER_PIOLLASTRESSNHFAST_API_H
#define _CODER_PIOLLASTRESSNHFAST_API_H

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
void PiollaStressNHFast(real_T F[9], real_T C100, real_T K, real_T S[9],
                        real_T D[36], real_T *P);

void PiollaStressNHFast_api(const mxArray *const prhs[3], int32_T nlhs,
                            const mxArray *plhs[3]);

void PiollaStressNHFast_atexit(void);

void PiollaStressNHFast_initialize(void);

void PiollaStressNHFast_terminate(void);

void PiollaStressNHFast_xil_shutdown(void);

void PiollaStressNHFast_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_PiollaStressNHFast_api.h) */
