/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlarfg.c
 *
 * Code generation for function 'xzlarfg'
 *
 */

/* Include files */
#include "xzlarfg.h"
#include "PolarDecompositionFast_rtwutil.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include <math.h>

/* Function Definitions */
double xzlarfg(int n, double *alpha1, double x[3])
{
  double beta1;
  double tau;
  double xnorm;
  int k;
  int knt;
  tau = 0.0;
  if (n > 0) {
    xnorm = xnrm2(n - 1, x);
    if (xnorm != 0.0) {
      beta1 = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        beta1 = -beta1;
      }
      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          for (k = 2; k <= n; k++) {
            x[k - 1] *= 9.9792015476736E+291;
          }
          beta1 *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(beta1) >= 1.0020841800044864E-292));
        beta1 = rt_hypotd_snf(*alpha1, xnrm2(n - 1, x));
        if (*alpha1 >= 0.0) {
          beta1 = -beta1;
        }
        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0 / (*alpha1 - beta1);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= xnorm;
        }
        for (k = 0; k <= knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }
        *alpha1 = beta1;
      } else {
        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0 / (*alpha1 - beta1);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= xnorm;
        }
        *alpha1 = beta1;
      }
    }
  }
  return tau;
}

/* End of code generation (xzlarfg.c) */
