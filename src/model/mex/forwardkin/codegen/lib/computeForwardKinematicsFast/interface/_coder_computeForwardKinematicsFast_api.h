/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_computeForwardKinematicsFast_api.h
 *
 * Code generation for function 'computeForwardKinematicsFast'
 *
 */

#ifndef _CODER_COMPUTEFORWARDKINEMATICSFAST_API_H
#define _CODER_COMPUTEFORWARDKINEMATICSFAST_API_H

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
void c_computeForwardKinematicsFast_(const mxArray *const prhs[7], int32_T nlhs,
                                     const mxArray *plhs[2]);

void computeForwardKinematicsFast(emxArray_real_T *x, real_T ds, real_T p0[3],
                                  real_T Phi0[9], emxArray_real_T *xia0,
                                  emxArray_real_T *Th, emxArray_real_T *Ba,
                                  emxArray_real_T *gtmp, emxArray_real_T *Jtmp);

void computeForwardKinematicsFast_atexit(void);

void computeForwardKinematicsFast_initialize(void);

void computeForwardKinematicsFast_terminate(void);

void computeForwardKinematicsFast_xil_shutdown(void);

void computeForwardKinematicsFast_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_computeForwardKinematicsFast_api.h) */
