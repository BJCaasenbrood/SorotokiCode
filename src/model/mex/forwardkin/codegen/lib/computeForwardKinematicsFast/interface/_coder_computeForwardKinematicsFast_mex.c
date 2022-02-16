/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_computeForwardKinematicsFast_mex.c
 *
 * Code generation for function 'computeForwardKinematicsFast'
 *
 */

/* Include files */
#include "_coder_computeForwardKinematicsFast_mex.h"
#include "_coder_computeForwardKinematicsFast_api.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&computeForwardKinematicsFast_atexit);
  /* Module initialization. */
  computeForwardKinematicsFast_initialize();
  /* Dispatch the entry-point. */
  unsafe_computeForwardKinematicsFast_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  computeForwardKinematicsFast_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_computeForwardKinematicsFast_mexFunction(int32_T nlhs,
                                                     mxArray *plhs[2],
                                                     int32_T nrhs,
                                                     const mxArray *prhs[7])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 7) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 7, 4,
                        28, "computeForwardKinematicsFast");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 28,
                        "computeForwardKinematicsFast");
  }
  /* Call the function. */
  c_computeForwardKinematicsFast_(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_computeForwardKinematicsFast_mex.c) */
