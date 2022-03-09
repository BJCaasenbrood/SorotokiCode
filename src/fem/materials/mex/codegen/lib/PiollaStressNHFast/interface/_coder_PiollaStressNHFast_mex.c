/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PiollaStressNHFast_mex.c
 *
 * Code generation for function 'PiollaStressNHFast'
 *
 */

/* Include files */
#include "_coder_PiollaStressNHFast_mex.h"
#include "_coder_PiollaStressNHFast_api.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&PiollaStressNHFast_atexit);
  /* Module initialization. */
  PiollaStressNHFast_initialize();
  /* Dispatch the entry-point. */
  unsafe_PiollaStressNHFast_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  PiollaStressNHFast_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_PiollaStressNHFast_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                           int32_T nrhs, const mxArray *prhs[3])
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
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4,
                        18, "PiollaStressNHFast");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 18,
                        "PiollaStressNHFast");
  }
  /* Call the function. */
  PiollaStressNHFast_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_PiollaStressNHFast_mex.c) */
