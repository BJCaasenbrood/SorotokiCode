/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PolarDecompositionFast.c
 *
 * Code generation for function 'PolarDecompositionFast'
 *
 */

/* Include files */
#include "PolarDecompositionFast.h"
#include "PolarDecompositionFast_data.h"
#include "eigHermitianStandard.h"
#include "rt_nonfinite.h"
#include "sqrt.h"
#include "xzggev.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void PolarDecompositionFast(const double F[9], creal_T R[9], creal_T S[9],
                            creal_T V[9])
{
  creal_T At[9];
  creal_T b_At[9];
  creal_T b_F[9];
  creal_T beta1[3];
  creal_T lambda[3];
  double C[9];
  double D[9];
  double b_V[9];
  double absxk;
  double colnorm;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double im;
  double re;
  double scale;
  double t;
  int coltop;
  int exitg1;
  int i;
  int i1;
  int ibmat;
  int k;
  int kend_tmp;
  boolean_T exitg2;
  boolean_T p;
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      C[i + 3 * i1] = (F[3 * i] * F[3 * i1] + F[3 * i + 1] * F[3 * i1 + 1]) +
                      F[3 * i + 2] * F[3 * i1 + 2];
    }
  }
  p = true;
  for (k = 0; k < 9; k++) {
    if ((!p) || (rtIsInf(C[k]) || rtIsNaN(C[k]))) {
      p = false;
    }
  }
  if (!p) {
    for (i = 0; i < 9; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
    }
    At[0].re = rtNaN;
    At[0].im = 0.0;
    At[4].re = rtNaN;
    At[4].im = 0.0;
    At[8].re = rtNaN;
    At[8].im = 0.0;
  } else {
    p = true;
    ibmat = 0;
    exitg2 = false;
    while ((!exitg2) && (ibmat < 3)) {
      k = 0;
      do {
        exitg1 = 0;
        if (k <= ibmat) {
          if (!(C[k + 3 * ibmat] == C[ibmat + 3 * k])) {
            p = false;
            exitg1 = 1;
          } else {
            k++;
          }
        } else {
          ibmat++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
    if (p) {
      eigHermitianStandard(C, b_V, D);
      for (i = 0; i < 9; i++) {
        V[i].re = b_V[i];
        V[i].im = 0.0;
        At[i].re = D[i];
        At[i].im = 0.0;
      }
    } else {
      for (i = 0; i < 9; i++) {
        At[i].re = C[i];
        At[i].im = 0.0;
      }
      xzggev(At, &k, lambda, beta1, V);
      for (coltop = 0; coltop <= 6; coltop += 3) {
        colnorm = 0.0;
        scale = 3.3121686421112381E-170;
        kend_tmp = coltop + 3;
        for (k = coltop + 1; k <= kend_tmp; k++) {
          absxk = fabs(V[k - 1].re);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }
          absxk = fabs(V[k - 1].im);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }
        }
        colnorm = scale * sqrt(colnorm);
        for (ibmat = coltop + 1; ibmat <= kend_tmp; ibmat++) {
          scale = V[ibmat - 1].re;
          absxk = V[ibmat - 1].im;
          if (absxk == 0.0) {
            re = scale / colnorm;
            im = 0.0;
          } else if (scale == 0.0) {
            re = 0.0;
            im = absxk / colnorm;
          } else {
            re = scale / colnorm;
            im = absxk / colnorm;
          }
          V[ibmat - 1].re = re;
          V[ibmat - 1].im = im;
        }
      }
      if (beta1[0].im == 0.0) {
        if (lambda[0].im == 0.0) {
          At[0].re = lambda[0].re / beta1[0].re;
          At[0].im = 0.0;
        } else if (lambda[0].re == 0.0) {
          At[0].re = 0.0;
          At[0].im = lambda[0].im / beta1[0].re;
        } else {
          At[0].re = lambda[0].re / beta1[0].re;
          At[0].im = lambda[0].im / beta1[0].re;
        }
      } else if (beta1[0].re == 0.0) {
        if (lambda[0].re == 0.0) {
          At[0].re = lambda[0].im / beta1[0].im;
          At[0].im = 0.0;
        } else if (lambda[0].im == 0.0) {
          At[0].re = 0.0;
          At[0].im = -(lambda[0].re / beta1[0].im);
        } else {
          At[0].re = lambda[0].im / beta1[0].im;
          At[0].im = -(lambda[0].re / beta1[0].im);
        }
      } else {
        colnorm = fabs(beta1[0].re);
        scale = fabs(beta1[0].im);
        if (colnorm > scale) {
          absxk = beta1[0].im / beta1[0].re;
          scale = beta1[0].re + absxk * beta1[0].im;
          At[0].re = (lambda[0].re + absxk * lambda[0].im) / scale;
          At[0].im = (lambda[0].im - absxk * lambda[0].re) / scale;
        } else if (scale == colnorm) {
          if (beta1[0].re > 0.0) {
            absxk = 0.5;
          } else {
            absxk = -0.5;
          }
          if (beta1[0].im > 0.0) {
            scale = 0.5;
          } else {
            scale = -0.5;
          }
          At[0].re = (lambda[0].re * absxk + lambda[0].im * scale) / colnorm;
          At[0].im = (lambda[0].im * absxk - lambda[0].re * scale) / colnorm;
        } else {
          absxk = beta1[0].re / beta1[0].im;
          scale = beta1[0].im + absxk * beta1[0].re;
          At[0].re = (absxk * lambda[0].re + lambda[0].im) / scale;
          At[0].im = (absxk * lambda[0].im - lambda[0].re) / scale;
        }
      }
      if (beta1[1].im == 0.0) {
        if (lambda[1].im == 0.0) {
          At[4].re = lambda[1].re / beta1[1].re;
          At[4].im = 0.0;
        } else if (lambda[1].re == 0.0) {
          At[4].re = 0.0;
          At[4].im = lambda[1].im / beta1[1].re;
        } else {
          At[4].re = lambda[1].re / beta1[1].re;
          At[4].im = lambda[1].im / beta1[1].re;
        }
      } else if (beta1[1].re == 0.0) {
        if (lambda[1].re == 0.0) {
          At[4].re = lambda[1].im / beta1[1].im;
          At[4].im = 0.0;
        } else if (lambda[1].im == 0.0) {
          At[4].re = 0.0;
          At[4].im = -(lambda[1].re / beta1[1].im);
        } else {
          At[4].re = lambda[1].im / beta1[1].im;
          At[4].im = -(lambda[1].re / beta1[1].im);
        }
      } else {
        colnorm = fabs(beta1[1].re);
        scale = fabs(beta1[1].im);
        if (colnorm > scale) {
          absxk = beta1[1].im / beta1[1].re;
          scale = beta1[1].re + absxk * beta1[1].im;
          At[4].re = (lambda[1].re + absxk * lambda[1].im) / scale;
          At[4].im = (lambda[1].im - absxk * lambda[1].re) / scale;
        } else if (scale == colnorm) {
          if (beta1[1].re > 0.0) {
            absxk = 0.5;
          } else {
            absxk = -0.5;
          }
          if (beta1[1].im > 0.0) {
            scale = 0.5;
          } else {
            scale = -0.5;
          }
          At[4].re = (lambda[1].re * absxk + lambda[1].im * scale) / colnorm;
          At[4].im = (lambda[1].im * absxk - lambda[1].re * scale) / colnorm;
        } else {
          absxk = beta1[1].re / beta1[1].im;
          scale = beta1[1].im + absxk * beta1[1].re;
          At[4].re = (absxk * lambda[1].re + lambda[1].im) / scale;
          At[4].im = (absxk * lambda[1].im - lambda[1].re) / scale;
        }
      }
      if (beta1[2].im == 0.0) {
        if (lambda[2].im == 0.0) {
          At[8].re = lambda[2].re / beta1[2].re;
          At[8].im = 0.0;
        } else if (lambda[2].re == 0.0) {
          At[8].re = 0.0;
          At[8].im = lambda[2].im / beta1[2].re;
        } else {
          At[8].re = lambda[2].re / beta1[2].re;
          At[8].im = lambda[2].im / beta1[2].re;
        }
      } else if (beta1[2].re == 0.0) {
        if (lambda[2].re == 0.0) {
          At[8].re = lambda[2].im / beta1[2].im;
          At[8].im = 0.0;
        } else if (lambda[2].im == 0.0) {
          At[8].re = 0.0;
          At[8].im = -(lambda[2].re / beta1[2].im);
        } else {
          At[8].re = lambda[2].im / beta1[2].im;
          At[8].im = -(lambda[2].re / beta1[2].im);
        }
      } else {
        colnorm = fabs(beta1[2].re);
        scale = fabs(beta1[2].im);
        if (colnorm > scale) {
          absxk = beta1[2].im / beta1[2].re;
          scale = beta1[2].re + absxk * beta1[2].im;
          At[8].re = (lambda[2].re + absxk * lambda[2].im) / scale;
          At[8].im = (lambda[2].im - absxk * lambda[2].re) / scale;
        } else if (scale == colnorm) {
          if (beta1[2].re > 0.0) {
            absxk = 0.5;
          } else {
            absxk = -0.5;
          }
          if (beta1[2].im > 0.0) {
            scale = 0.5;
          } else {
            scale = -0.5;
          }
          At[8].re = (lambda[2].re * absxk + lambda[2].im * scale) / colnorm;
          At[8].im = (lambda[2].im * absxk - lambda[2].re * scale) / colnorm;
        } else {
          absxk = beta1[2].re / beta1[2].im;
          scale = beta1[2].im + absxk * beta1[2].re;
          At[8].re = (absxk * lambda[2].re + lambda[2].im) / scale;
          At[8].im = (absxk * lambda[2].im - lambda[2].re) / scale;
        }
      }
    }
  }
  lambda[0] = At[0];
  lambda[1] = At[4];
  lambda[2] = At[8];
  b_sqrt(&lambda[0]);
  b_sqrt(&lambda[1]);
  b_sqrt(&lambda[2]);
  for (k = 0; k < 3; k++) {
    t = lambda[k].re;
    im = -lambda[k].im;
    if (im == 0.0) {
      re = 1.0 / t;
      im = 0.0;
    } else if (t == 0.0) {
      re = 0.0;
      im = -(1.0 / im);
    } else {
      colnorm = fabs(t);
      scale = fabs(im);
      if (colnorm > scale) {
        absxk = im / t;
        scale = t + absxk * im;
        re = (absxk * 0.0 + 1.0) / scale;
        im = (0.0 - absxk) / scale;
      } else if (scale == colnorm) {
        if (t > 0.0) {
          absxk = 0.5;
        } else {
          absxk = -0.5;
        }
        if (im > 0.0) {
          scale = 0.5;
        } else {
          scale = -0.5;
        }
        re = (absxk + 0.0 * scale) / colnorm;
        im = (0.0 * absxk - scale) / colnorm;
      } else {
        absxk = t / im;
        scale = im + absxk * t;
        re = absxk / scale;
        im = (absxk * 0.0 - 1.0) / scale;
      }
    }
    lambda[k].re = re;
    lambda[k].im = im;
    ibmat = k * 3;
    At[ibmat] = lambda[k];
    At[ibmat + 1] = lambda[k];
    At[ibmat + 2] = lambda[k];
  }
  for (i = 0; i < 9; i++) {
    d = At[i].re;
    d1 = At[i].im;
    d2 = V[i].im;
    d3 = V[i].re;
    im = d * d2 + d1 * d3;
    d = d * d3 - d1 * d2;
    At[i].re = d;
    At[i].im = im;
    b_F[i].re = F[i];
    b_F[i].im = 0.0;
  }
  for (i = 0; i < 3; i++) {
    d = At[i].re;
    d1 = At[i].im;
    d2 = At[i + 3].re;
    d3 = At[i + 3].im;
    d4 = At[i + 6].re;
    d5 = At[i + 6].im;
    for (i1 = 0; i1 < 3; i1++) {
      scale = V[i1].re;
      absxk = -V[i1].im;
      re = d * scale - d1 * absxk;
      im = d * absxk + d1 * scale;
      scale = V[i1 + 3].re;
      absxk = -V[i1 + 3].im;
      re += d2 * scale - d3 * absxk;
      im += d2 * absxk + d3 * scale;
      scale = V[i1 + 6].re;
      absxk = -V[i1 + 6].im;
      re += d4 * scale - d5 * absxk;
      im += d4 * absxk + d5 * scale;
      kend_tmp = i + 3 * i1;
      b_At[kend_tmp].re = re;
      b_At[kend_tmp].im = im;
    }
  }
  for (i = 0; i < 3; i++) {
    d = b_F[i].re;
    d1 = b_F[i].im;
    d2 = b_F[i + 3].re;
    d3 = b_F[i + 3].im;
    d4 = b_F[i + 6].re;
    d5 = b_F[i + 6].im;
    for (i1 = 0; i1 < 3; i1++) {
      ibmat = 3 * i1 + 1;
      k = 3 * i1 + 2;
      kend_tmp = i + 3 * i1;
      scale = b_At[3 * i1].re;
      absxk = b_At[3 * i1].im;
      t = b_At[ibmat].re;
      im = b_At[ibmat].im;
      colnorm = b_At[k].re;
      re = b_At[k].im;
      R[kend_tmp].re = ((d * scale - d1 * absxk) + (d2 * t - d3 * im)) +
                       (d4 * colnorm - d5 * re);
      R[kend_tmp].im = ((d * absxk + d1 * scale) + (d2 * im + d3 * t)) +
                       (d4 * re + d5 * colnorm);
    }
  }
  for (i = 0; i < 9; i++) {
    b_F[i].re = F[i];
    b_F[i].im = 0.0;
  }
  for (i = 0; i < 3; i++) {
    scale = R[3 * i].re;
    d = R[3 * i].im;
    k = 3 * i + 1;
    coltop = 3 * i + 2;
    for (i1 = 0; i1 < 3; i1++) {
      d1 = b_F[3 * i1].re;
      d2 = b_F[3 * i1].im;
      absxk = R[k].re;
      t = -R[k].im;
      ibmat = 3 * i1 + 1;
      d3 = b_F[ibmat].re;
      d4 = b_F[ibmat].im;
      re = (scale * d1 - -d * d2) + (absxk * d3 - t * d4);
      im = (scale * d2 + -d * d1) + (absxk * d4 + t * d3);
      absxk = R[coltop].re;
      t = -R[coltop].im;
      ibmat = 3 * i1 + 2;
      d1 = b_F[ibmat].re;
      d2 = b_F[ibmat].im;
      re += absxk * d1 - t * d2;
      im += absxk * d2 + t * d1;
      kend_tmp = i + 3 * i1;
      S[kend_tmp].re = re;
      S[kend_tmp].im = im;
    }
  }
  for (i = 0; i < 9; i++) {
    b_F[i].re = F[i];
    b_F[i].im = 0.0;
  }
  for (i = 0; i < 3; i++) {
    d = b_F[i].re;
    d1 = b_F[i].im;
    d2 = b_F[i + 3].re;
    d3 = b_F[i + 3].im;
    d4 = b_F[i + 6].re;
    d5 = b_F[i + 6].im;
    for (i1 = 0; i1 < 3; i1++) {
      absxk = R[i1].re;
      t = -R[i1].im;
      re = d * absxk - d1 * t;
      im = d * t + d1 * absxk;
      absxk = R[i1 + 3].re;
      t = -R[i1 + 3].im;
      re += d2 * absxk - d3 * t;
      im += d2 * t + d3 * absxk;
      absxk = R[i1 + 6].re;
      t = -R[i1 + 6].im;
      re += d4 * absxk - d5 * t;
      im += d4 * t + d5 * absxk;
      kend_tmp = i + 3 * i1;
      V[kend_tmp].re = re;
      V[kend_tmp].im = im;
    }
  }
}

/* End of code generation (PolarDecompositionFast.c) */
