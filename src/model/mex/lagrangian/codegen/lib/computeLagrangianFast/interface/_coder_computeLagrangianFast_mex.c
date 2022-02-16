/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_computeLagrangianFast_mex.c
 *
 * Code generation for function 'computeLagrangianFast'
 *
 */

/* Include files */
#include "_coder_computeLagrangianFast_mex.h"
#include "_coder_computeLagrangianFast_api.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&computeLagrangianFast_atexit);
  /* Module initialization. */
  computeLagrangianFast_initialize();
  /* Dispatch the entry-point. */
  unsafe_computeLagrangianFast_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  computeLagrangianFast_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_computeLagrangianFast_mexFunction(int32_T nlhs, mxArray *plhs[10],
                                              int32_T nrhs,
                                              const mxArray *prhs[11])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[10];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 11) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 11, 4,
                        21, "computeLagrangianFast");
  }
  if (nlhs > 10) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 21,
                        "computeLagrangianFast");
  }
  /* Call the function. */
  computeLagrangianFast_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_computeLagrangianFast_mex.c) */
