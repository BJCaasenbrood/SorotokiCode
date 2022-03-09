/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xztgevc.c
 *
 * Code generation for function 'xztgevc'
 *
 */

/* Include files */
#include "xztgevc.h"
#include "PolarDecompositionFast_data.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void xztgevc(const creal_T A[9], creal_T V[9])
{
  creal_T work1[3];
  creal_T work2[3];
  double rworka[3];
  double acoeff;
  double ai;
  double anorm;
  double ascale;
  double brm;
  double d_im;
  double d_re;
  double dmin;
  double salpha_im;
  double salpha_re;
  double scale;
  double temp;
  double xmx;
  double y;
  int b_i;
  int b_j;
  int i;
  int j;
  int je;
  int re_tmp;
  int x_tmp_tmp_tmp;
  boolean_T lscalea;
  boolean_T lscaleb;
  rworka[0] = 0.0;
  rworka[1] = 0.0;
  rworka[2] = 0.0;
  anorm = fabs(A[0].re) + fabs(A[0].im);
  for (j = 0; j < 2; j++) {
    for (i = 0; i <= j; i++) {
      re_tmp = i + 3 * (j + 1);
      rworka[j + 1] += fabs(A[re_tmp].re) + fabs(A[re_tmp].im);
    }
    re_tmp = (j + 3 * (j + 1)) + 1;
    y = rworka[j + 1] + (fabs(A[re_tmp].re) + fabs(A[re_tmp].im));
    if (y > anorm) {
      anorm = y;
    }
  }
  y = anorm;
  if (2.2250738585072014E-308 > anorm) {
    y = 2.2250738585072014E-308;
  }
  ascale = 1.0 / y;
  for (je = 0; je < 3; je++) {
    x_tmp_tmp_tmp = 3 * (2 - je);
    re_tmp = (x_tmp_tmp_tmp - je) + 2;
    xmx = A[re_tmp].re;
    scale = A[re_tmp].im;
    y = (fabs(xmx) + fabs(scale)) * ascale;
    if (1.0 > y) {
      y = 1.0;
    }
    temp = 1.0 / y;
    salpha_re = ascale * (temp * xmx);
    salpha_im = ascale * (temp * scale);
    acoeff = temp * ascale;
    if ((temp >= 2.2250738585072014E-308) &&
        (acoeff < 3.0062525400134592E-292)) {
      lscalea = true;
    } else {
      lscalea = false;
    }
    xmx = fabs(salpha_re) + fabs(salpha_im);
    if ((xmx >= 2.2250738585072014E-308) && (xmx < 3.0062525400134592E-292)) {
      lscaleb = true;
    } else {
      lscaleb = false;
    }
    scale = 1.0;
    if (lscalea) {
      y = anorm;
      if (3.3264005158911995E+291 < anorm) {
        y = 3.3264005158911995E+291;
      }
      scale = 3.0062525400134592E-292 / temp * y;
    }
    if (lscaleb) {
      y = 3.0062525400134592E-292 / xmx;
      if (y > scale) {
        scale = y;
      }
    }
    if (lscalea || lscaleb) {
      y = acoeff;
      if (1.0 > acoeff) {
        y = 1.0;
      }
      if (xmx > y) {
        y = xmx;
      }
      y = 1.0 / (2.2250738585072014E-308 * y);
      if (y < scale) {
        scale = y;
      }
      if (lscalea) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }
      salpha_re *= scale;
      salpha_im *= scale;
    }
    work1[0].re = 0.0;
    work1[0].im = 0.0;
    work1[1].re = 0.0;
    work1[1].im = 0.0;
    work1[2].re = 0.0;
    work1[2].im = 0.0;
    work1[2 - je].re = 1.0;
    work1[2 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoeff * anorm;
    y = 2.2204460492503131E-16 * (fabs(salpha_re) + fabs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }
    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }
    b_i = 1 - je;
    for (i = 0; i <= b_i; i++) {
      re_tmp = i + x_tmp_tmp_tmp;
      work1[i].re = acoeff * A[re_tmp].re;
      work1[i].im = acoeff * A[re_tmp].im;
    }
    work1[2 - je].re = 1.0;
    work1[2 - je].im = 0.0;
    b_i = (int)(((-1.0 - ((-(double)je + 3.0) - 1.0)) + 1.0) / -1.0);
    for (j = 0; j < b_i; j++) {
      b_j = 1 - (je + j);
      re_tmp = b_j + 3 * b_j;
      d_re = acoeff * A[re_tmp].re - salpha_re;
      d_im = acoeff * A[re_tmp].im - salpha_im;
      if (fabs(d_re) + fabs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0;
      }
      brm = fabs(d_re);
      y = fabs(d_im);
      xmx = brm + y;
      if (xmx < 1.0) {
        scale = fabs(work1[b_j].re) + fabs(work1[b_j].im);
        if (scale >= 1.4980776123852632E+307 * xmx) {
          temp = 1.0 / scale;
          re_tmp = 2 - je;
          for (i = 0; i <= re_tmp; i++) {
            work1[i].re *= temp;
            work1[i].im *= temp;
          }
        }
      }
      temp = -work1[b_j].re;
      ai = -work1[b_j].im;
      if (d_im == 0.0) {
        if (ai == 0.0) {
          y = temp / d_re;
          xmx = 0.0;
        } else if (temp == 0.0) {
          y = 0.0;
          xmx = ai / d_re;
        } else {
          y = temp / d_re;
          xmx = ai / d_re;
        }
      } else if (d_re == 0.0) {
        if (temp == 0.0) {
          y = ai / d_im;
          xmx = 0.0;
        } else if (ai == 0.0) {
          y = 0.0;
          xmx = -(temp / d_im);
        } else {
          y = ai / d_im;
          xmx = -(temp / d_im);
        }
      } else if (brm > y) {
        scale = d_im / d_re;
        xmx = d_re + scale * d_im;
        y = (temp + scale * ai) / xmx;
        xmx = (ai - scale * temp) / xmx;
      } else if (y == brm) {
        if (d_re > 0.0) {
          scale = 0.5;
        } else {
          scale = -0.5;
        }
        if (d_im > 0.0) {
          xmx = 0.5;
        } else {
          xmx = -0.5;
        }
        y = (temp * scale + ai * xmx) / brm;
        xmx = (ai * scale - temp * xmx) / brm;
      } else {
        scale = d_re / d_im;
        xmx = d_im + scale * d_re;
        y = (scale * temp + ai) / xmx;
        xmx = (scale * ai - temp) / xmx;
      }
      work1[b_j].re = y;
      work1[b_j].im = xmx;
      if (b_j + 1 > 1) {
        xmx = fabs(work1[1].re) + fabs(work1[1].im);
        if (xmx > 1.0) {
          temp = 1.0 / xmx;
          if (acoeff * rworka[1] >= 1.4980776123852632E+307 * temp) {
            re_tmp = 2 - je;
            for (i = 0; i <= re_tmp; i++) {
              work1[i].re *= temp;
              work1[i].im *= temp;
            }
          }
        }
        d_re = acoeff * work1[1].re;
        d_im = acoeff * work1[1].im;
        work1[0].re += d_re * A[3].re - d_im * A[3].im;
        work1[0].im += d_re * A[3].im + d_im * A[3].re;
      }
    }
    work2[0].re = 0.0;
    work2[0].im = 0.0;
    work2[1].re = 0.0;
    work2[1].im = 0.0;
    work2[2].re = 0.0;
    work2[2].im = 0.0;
    b_i = 2 - je;
    for (i = 0; i <= b_i; i++) {
      xmx = work1[i].re;
      scale = work1[i].im;
      work2[0].re += V[3 * i].re * xmx - V[3 * i].im * scale;
      work2[0].im += V[3 * i].re * scale + V[3 * i].im * xmx;
      re_tmp = 3 * i + 1;
      work2[1].re += V[re_tmp].re * xmx - V[re_tmp].im * scale;
      work2[1].im += V[re_tmp].re * scale + V[re_tmp].im * xmx;
      re_tmp = 3 * i + 2;
      work2[2].re += V[re_tmp].re * xmx - V[re_tmp].im * scale;
      work2[2].im += V[re_tmp].re * scale + V[re_tmp].im * xmx;
    }
    xmx = fabs(work2[0].re) + fabs(work2[0].im);
    y = fabs(work2[1].re) + fabs(work2[1].im);
    if (y > xmx) {
      xmx = y;
    }
    y = fabs(work2[2].re) + fabs(work2[2].im);
    if (y > xmx) {
      xmx = y;
    }
    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      V[x_tmp_tmp_tmp].re = temp * work2[0].re;
      V[x_tmp_tmp_tmp].im = temp * work2[0].im;
      V[x_tmp_tmp_tmp + 1].re = temp * work2[1].re;
      V[x_tmp_tmp_tmp + 1].im = temp * work2[1].im;
      V[x_tmp_tmp_tmp + 2].re = temp * work2[2].re;
      V[x_tmp_tmp_tmp + 2].im = temp * work2[2].im;
    } else {
      V[x_tmp_tmp_tmp].re = 0.0;
      V[x_tmp_tmp_tmp].im = 0.0;
      V[x_tmp_tmp_tmp + 1].re = 0.0;
      V[x_tmp_tmp_tmp + 1].im = 0.0;
      V[x_tmp_tmp_tmp + 2].re = 0.0;
      V[x_tmp_tmp_tmp + 2].im = 0.0;
    }
  }
}

/* End of code generation (xztgevc.c) */
