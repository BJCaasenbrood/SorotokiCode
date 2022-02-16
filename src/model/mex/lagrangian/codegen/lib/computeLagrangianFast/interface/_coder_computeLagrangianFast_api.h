/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_computeLagrangianFast_api.h
 *
 * Code generation for function 'computeLagrangianFast'
 *
 */

#ifndef _CODER_COMPUTELAGRANGIANFAST_API_H
#define _CODER_COMPUTELAGRANGIANFAST_API_H

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
void computeLagrangianFast(emxArray_real_T *x, emxArray_real_T *dx, real_T ds,
                           real_T p0[3], real_T Phi0[9], emxArray_real_T *xia0,
                           emxArray_real_T *Th, emxArray_real_T *Ba,
                           real_T Ktt[36], real_T Mtt[36], real_T Zeta,
                           emxArray_real_T *M, emxArray_real_T *C,
                           emxArray_real_T *K, emxArray_real_T *R,
                           emxArray_real_T *G, real_T p[3], real_T Phi[9],
                           emxArray_real_T *J, real_T *Vg, real_T *Kin);

void computeLagrangianFast_api(const mxArray *const prhs[11], int32_T nlhs,
                               const mxArray *plhs[10]);

void computeLagrangianFast_atexit(void);

void computeLagrangianFast_initialize(void);

void computeLagrangianFast_terminate(void);

void computeLagrangianFast_xil_shutdown(void);

void computeLagrangianFast_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_computeLagrangianFast_api.h) */
