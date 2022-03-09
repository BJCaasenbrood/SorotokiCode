/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_NonlinearStrainOperatorFast_api.h
 *
 * Code generation for function 'NonlinearStrainOperatorFast'
 *
 */

#ifndef _CODER_NONLINEARSTRAINOPERATORFAST_API_H
#define _CODER_NONLINEARSTRAINOPERATORFAST_API_H

/* Include files */
#include "emlrt.h"
#include "tmwtypes.h"
#include <string.h>

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void NonlinearStrainOperatorFast(emxArray_real_T *N, emxArray_real_T *dNdx,
                                 real_T F[9], emxArray_real_T *Bn,
                                 emxArray_real_T *Bg, emxArray_real_T *NN,
                                 real_T *tau);

void NonlinearStrainOperatorFast_api(const mxArray *const prhs[3], int32_T nlhs,
                                     const mxArray *plhs[4]);

void NonlinearStrainOperatorFast_atexit(void);

void NonlinearStrainOperatorFast_initialize(void);

void NonlinearStrainOperatorFast_terminate(void);

void NonlinearStrainOperatorFast_xil_shutdown(void);

void NonlinearStrainOperatorFast_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_NonlinearStrainOperatorFast_api.h) */
