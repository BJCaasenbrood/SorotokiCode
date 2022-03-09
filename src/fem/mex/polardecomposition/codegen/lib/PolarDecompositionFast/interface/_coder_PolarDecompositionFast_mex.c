/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PolarDecompositionFast_mex.c
 *
 * Code generation for function 'PolarDecompositionFast'
 *
 */

/* Include files */
#include "_coder_PolarDecompositionFast_mex.h"
#include "_coder_PolarDecompositionFast_api.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&PolarDecompositionFast_atexit);
  /* Module initialization. */
  PolarDecompositionFast_initialize();
  /* Dispatch the entry-point. */
  unsafe_PolarDecompositionFast_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  PolarDecompositionFast_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_PolarDecompositionFast_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                               int32_T nrhs,
                                               const mxArray *prhs[1])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[3];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 1, 4,
                        22, "PolarDecompositionFast");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "PolarDecompositionFast");
  }
  /* Call the function. */
  PolarDecompositionFast_api(prhs[0], nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_PolarDecompositionFast_mex.c) */
