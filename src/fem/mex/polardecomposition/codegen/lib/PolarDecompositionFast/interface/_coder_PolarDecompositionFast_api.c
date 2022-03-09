/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PolarDecompositionFast_api.c
 *
 * Code generation for function 'PolarDecompositionFast'
 *
 */

/* Include files */
#include "_coder_PolarDecompositionFast_api.h"
#include "_coder_PolarDecompositionFast_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131610U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "PolarDecompositionFast",                             /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[9];

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[9];

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *F,
                                 const char_T *identifier))[9];

static const mxArray *emlrt_marshallOut(const emlrtStack *sp,
                                        const creal_T u[9]);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[9]
{
  real_T(*y)[9];
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[9]
{
  static const int32_T dims[2] = {3, 3};
  real_T(*ret)[9];
  emlrtCheckBuiltInR2012b((emlrtCTX)sp, msgId, src, (const char_T *)"double",
                          false, 2U, (void *)&dims[0]);
  ret = (real_T(*)[9])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *F,
                                 const char_T *identifier))[9]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[9];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(F), &thisId);
  emlrtDestroyArray(&F);
  return y;
}

static const mxArray *emlrt_marshallOut(const emlrtStack *sp,
                                        const creal_T u[9])
{
  static const int32_T iv[2] = {3, 3};
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS,
                              mxCOMPLEX);
  emlrtExportNumericArrayR2013b((emlrtCTX)sp, m, (void *)&u[0], 8);
  emlrtAssign(&y, m);
  return y;
}

void PolarDecompositionFast_api(const mxArray *prhs, int32_T nlhs,
                                const mxArray *plhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  creal_T R[9];
  creal_T S[9];
  creal_T V[9];
  real_T(*F)[9];
  st.tls = emlrtRootTLSGlobal;
  /* Marshall function inputs */
  F = emlrt_marshallIn(&st, emlrtAlias(prhs), "F");
  /* Invoke the target function */
  PolarDecompositionFast(*F, R, S, V);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&st, R);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(&st, S);
  }
  if (nlhs > 2) {
    plhs[2] = emlrt_marshallOut(&st, V);
  }
}

void PolarDecompositionFast_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  PolarDecompositionFast_xil_terminate();
  PolarDecompositionFast_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void PolarDecompositionFast_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void PolarDecompositionFast_terminate(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_PolarDecompositionFast_api.c) */
