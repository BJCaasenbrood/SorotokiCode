/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eigHermitianStandard.c
 *
 * Code generation for function 'eigHermitianStandard'
 *
 */

/* Include files */
#include "eigHermitianStandard.h"
#include "PolarDecompositionFast_rtwutil.h"
#include "rt_nonfinite.h"
#include "xgerc.h"
#include "xhseqr.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static int div_nde_s32_floor(int numerator, int denominator);

/* Function Definitions */
static int div_nde_s32_floor(int numerator, int denominator)
{
  int b_numerator;
  if (((numerator < 0) != (denominator < 0)) &&
      (numerator % denominator != 0)) {
    b_numerator = -1;
  } else {
    b_numerator = 0;
  }
  return numerator / denominator + b_numerator;
}

void eigHermitianStandard(const double A[9], double V[9], double D[9])
{
  double work[3];
  double tau[2];
  double b_alpha1_tmp;
  double beta1;
  double xnorm;
  int alpha1_tmp;
  int b_i;
  int c_i;
  int exitg1;
  int i;
  int i1;
  int ia;
  int in;
  int iv0_tmp;
  int ix0;
  int k;
  int knt;
  int lastc;
  int lastv;
  boolean_T exitg2;
  boolean_T p;
  p = true;
  for (k = 0; k < 9; k++) {
    if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
      p = false;
    }
  }
  if (!p) {
    for (i = 0; i < 9; i++) {
      V[i] = rtNaN;
    }
    b_i = 2;
    for (k = 0; k < 2; k++) {
      if (b_i <= 3) {
        memset(&V[(k * 3 + b_i) + -1], 0, (4 - b_i) * sizeof(double));
      }
      b_i++;
    }
    for (i = 0; i < 9; i++) {
      D[i] = rtNaN;
    }
  } else {
    memcpy(&D[0], &A[0], 9U * sizeof(double));
    work[0] = 0.0;
    work[1] = 0.0;
    work[2] = 0.0;
    for (c_i = 0; c_i < 2; c_i++) {
      b_i = c_i * 3 + 2;
      in = (c_i + 1) * 3;
      alpha1_tmp = (c_i + 3 * c_i) + 1;
      b_alpha1_tmp = D[alpha1_tmp];
      ix0 = b_i + 1;
      tau[c_i] = 0.0;
      xnorm = 0.0;
      if (1 - c_i >= 1) {
        xnorm = fabs(D[b_i]);
      }
      if (xnorm != 0.0) {
        beta1 = rt_hypotd_snf(b_alpha1_tmp, xnorm);
        if (b_alpha1_tmp >= 0.0) {
          beta1 = -beta1;
        }
        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i = (b_i - c_i) + 1;
          do {
            knt++;
            for (k = ix0; k <= i; k++) {
              D[k - 1] *= 9.9792015476736E+291;
            }
            beta1 *= 9.9792015476736E+291;
            b_alpha1_tmp *= 9.9792015476736E+291;
          } while (!(fabs(beta1) >= 1.0020841800044864E-292));
          xnorm = 0.0;
          if (1 - c_i >= 1) {
            xnorm = fabs(D[b_i]);
          }
          beta1 = rt_hypotd_snf(b_alpha1_tmp, xnorm);
          if (b_alpha1_tmp >= 0.0) {
            beta1 = -beta1;
          }
          tau[c_i] = (beta1 - b_alpha1_tmp) / beta1;
          xnorm = 1.0 / (b_alpha1_tmp - beta1);
          for (k = ix0; k <= i; k++) {
            D[k - 1] *= xnorm;
          }
          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }
          b_alpha1_tmp = beta1;
        } else {
          tau[c_i] = (beta1 - b_alpha1_tmp) / beta1;
          xnorm = 1.0 / (b_alpha1_tmp - beta1);
          i = (b_i - c_i) + 1;
          for (k = ix0; k <= i; k++) {
            D[k - 1] *= xnorm;
          }
          b_alpha1_tmp = beta1;
        }
      }
      D[alpha1_tmp] = 1.0;
      iv0_tmp = c_i + b_i;
      ix0 = in + 1;
      if (tau[c_i] != 0.0) {
        lastv = 1 - c_i;
        b_i = iv0_tmp - c_i;
        while ((lastv + 1 > 0) && (D[b_i] == 0.0)) {
          lastv--;
          b_i--;
        }
        lastc = 3;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = in + lastc;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= knt + lastv * 3) {
              if (D[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia += 3;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);
          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = -1;
        lastc = 0;
      }
      if (lastv + 1 > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            memset(&work[0], 0, lastc * sizeof(double));
          }
          knt = iv0_tmp - 1;
          i = (in + 3 * lastv) + 1;
          for (b_i = ix0; b_i <= i; b_i += 3) {
            i1 = (b_i + lastc) - 1;
            for (ia = b_i; ia <= i1; ia++) {
              k = ia - b_i;
              work[k] += D[ia - 1] * D[knt];
            }
            knt++;
          }
        }
        if (!(-tau[c_i] == 0.0)) {
          knt = in;
          for (k = 0; k <= lastv; k++) {
            xnorm = D[(iv0_tmp + k) - 1];
            if (xnorm != 0.0) {
              xnorm *= -tau[c_i];
              i = knt + 1;
              i1 = lastc + knt;
              for (b_i = i; b_i <= i1; b_i++) {
                D[b_i - 1] += work[(b_i - knt) - 1] * xnorm;
              }
            }
            knt += 3;
          }
        }
      }
      ix0 = (c_i + in) + 2;
      if (tau[c_i] != 0.0) {
        lastv = 2 - c_i;
        b_i = iv0_tmp - c_i;
        while ((lastv > 0) && (D[b_i] == 0.0)) {
          lastv--;
          b_i--;
        }
        lastc = 2 - c_i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = ix0 + (lastc - 1) * 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (D[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);
          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }
      if (lastv > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            memset(&work[0], 0, lastc * sizeof(double));
          }
          i = ix0 + 3 * (lastc - 1);
          for (b_i = ix0; b_i <= i; b_i += 3) {
            xnorm = 0.0;
            i1 = (b_i + lastv) - 1;
            for (ia = b_i; ia <= i1; ia++) {
              xnorm += D[ia - 1] * D[((iv0_tmp + ia) - b_i) - 1];
            }
            k = div_nde_s32_floor(b_i - ix0, 3);
            work[k] += xnorm;
          }
        }
        xgerc(lastv, lastc, -tau[c_i], iv0_tmp, work, D, ix0);
      }
      D[alpha1_tmp] = b_alpha1_tmp;
    }
    memcpy(&V[0], &D[0], 9U * sizeof(double));
    for (k = 1; k >= 0; k--) {
      ia = (k + 1) * 3;
      for (c_i = 0; c_i <= k; c_i++) {
        V[ia + c_i] = 0.0;
      }
      i = k + 3;
      for (c_i = i; c_i < 4; c_i++) {
        V[ia + 2] = V[ia - 1];
      }
    }
    V[1] = 0.0;
    V[2] = 0.0;
    V[0] = 1.0;
    work[0] = 0.0;
    work[1] = 0.0;
    work[2] = 0.0;
    for (c_i = 1; c_i >= 0; c_i--) {
      knt = (c_i + c_i * 3) + 8;
      if (c_i + 1 < 2) {
        V[knt - 4] = 1.0;
        if (tau[c_i] != 0.0) {
          lastv = 2;
          b_i = knt;
          while ((lastv > 0) && (V[b_i - 3] == 0.0)) {
            lastv--;
            b_i--;
          }
          lastc = 1;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (V[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc = 0;
              exitg1 = 1;
            }
          } while (exitg1 == 0);
        } else {
          lastv = 0;
          lastc = 0;
        }
        if (lastv > 0) {
          if (lastc != 0) {
            work[0] = 0.0;
            for (b_i = knt; b_i <= knt; b_i += 3) {
              xnorm = 0.0;
              i = (b_i + lastv) - 1;
              for (ia = b_i; ia <= i; ia++) {
                xnorm += V[ia - 1] * V[((knt + ia) - b_i) - 4];
              }
              k = div_nde_s32_floor(b_i - knt, 3);
              work[k] += xnorm;
            }
          }
          xgerc(lastv, lastc, -tau[c_i], knt - 3, work, V, knt);
        }
        ix0 = knt - 2;
        for (k = ix0; k <= ix0; k++) {
          V[k - 1] *= -tau[c_i];
        }
      }
      V[knt - 4] = 1.0 - tau[c_i];
      if (0 <= c_i - 1) {
        V[knt - 5] = 0.0;
      }
    }
    xhseqr(D, V);
  }
  for (k = 0; k < 2; k++) {
    D[(k + 3 * k) + 1] = 0.0;
    for (c_i = 0; c_i <= k; c_i++) {
      D[c_i + 3 * (k + 1)] = 0.0;
    }
  }
}

/* End of code generation (eigHermitianStandard.c) */
