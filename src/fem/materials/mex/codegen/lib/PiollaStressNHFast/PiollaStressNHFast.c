/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PiollaStressNHFast.c
 *
 * Code generation for function 'PiollaStressNHFast'
 *
 */

/* Include files */
#include "PiollaStressNHFast.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static double rt_powd_snf(double u0, double u1)
{
  double d;
  double d1;
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

void PiollaStressNHFast(const double F[9], double C100, double K, double S[9],
                        double D[36], double *P)
{
  static const signed char b[6] = {2, 2, 2, 0, 0, 0};
  double I3EE[36];
  double J1EE_tmp[36];
  double b_I3E[36];
  double C[9];
  double x[9];
  double I3E[6];
  double J1E[6];
  double Se[6];
  double I1;
  double a_tmp;
  double a_tmp_tmp;
  double b_a_tmp;
  double c_a_tmp;
  double d;
  double d1;
  double d_a_tmp;
  double e_a_tmp;
  double s;
  double smax;
  int a;
  int b_tmp;
  int i;
  int j;
  int jA;
  int jp1j;
  int k;
  int mmj_tmp;
  signed char ipiv[3];
  boolean_T isodd;
  /* Se = 2nd PK stress [S11, S22, S33, S12, S23, S13]; */
  /*  if (nargin > 2) && Robustness */
  /*  %Nu0 = Nu0*0.75; */
  /*  end */
  /* C100 = NeoHookeanMaterial.C10;  */
  /* K    = NeoHookeanMaterial.D1; */
  for (i = 0; i < 3; i++) {
    for (jp1j = 0; jp1j < 3; jp1j++) {
      C[i + 3 * jp1j] =
          (F[3 * i] * F[3 * jp1j] + F[3 * i + 1] * F[3 * jp1j + 1]) +
          F[3 * i + 2] * F[3 * jp1j + 2];
    }
  }
  I1 = (C[0] + C[4]) + C[8];
  memcpy(&x[0], &C[0], 9U * sizeof(double));
  ipiv[0] = 1;
  ipiv[1] = 2;
  for (j = 0; j < 2; j++) {
    mmj_tmp = 1 - j;
    b_tmp = j << 2;
    jp1j = b_tmp + 2;
    jA = 3 - j;
    a = 0;
    smax = fabs(x[b_tmp]);
    for (k = 2; k <= jA; k++) {
      s = fabs(x[(b_tmp + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (x[b_tmp + a] != 0.0) {
      if (a != 0) {
        jA = j + a;
        ipiv[j] = (signed char)(jA + 1);
        smax = x[j];
        x[j] = x[jA];
        x[jA] = smax;
        smax = x[j + 3];
        x[j + 3] = x[jA + 3];
        x[jA + 3] = smax;
        smax = x[j + 6];
        x[j + 6] = x[jA + 6];
        x[jA + 6] = smax;
      }
      i = (b_tmp - j) + 3;
      for (jA = jp1j; jA <= i; jA++) {
        x[jA - 1] /= x[b_tmp];
      }
    }
    jA = b_tmp;
    for (a = 0; a <= mmj_tmp; a++) {
      smax = x[(b_tmp + a * 3) + 3];
      if (smax != 0.0) {
        i = jA + 5;
        jp1j = (jA - j) + 6;
        for (k = i; k <= jp1j; k++) {
          x[k - 1] += x[((b_tmp + k) - jA) - 4] * -smax;
        }
      }
      jA += 3;
    }
  }
  isodd = (ipiv[0] > 1);
  smax = x[0] * x[4] * x[8];
  if (ipiv[1] > 2) {
    isodd = !isodd;
  }
  if (isodd) {
    smax = -smax;
  }
  s = sqrt(smax);
  /*  */
  I3E[0] = 2.0 * (C[4] * C[8] - C[7] * C[7]);
  I3E[1] = 2.0 * (C[0] * C[8] - C[6] * C[6]);
  I3E[2] = 2.0 * (C[0] * C[4] - C[3] * C[3]);
  I3E[3] = 2.0 * (C[6] * C[7] - C[3] * C[8]);
  I3E[4] = 2.0 * (C[3] * C[6] - C[0] * C[7]);
  I3E[5] = 2.0 * (C[3] * C[7] - C[4] * C[6]);
  /*  */
  a_tmp = rt_powd_snf(smax, -0.33333333333333331);
  b_a_tmp = rt_powd_snf(smax, -1.3333333333333333);
  c_a_tmp = 0.33333333333333331 * I1 * b_a_tmp;
  a_tmp_tmp = rt_powd_snf(smax, -0.5);
  d_a_tmp = 0.5 * a_tmp_tmp;
  /*  */
  /*  */
  *P = C100 * (I1 * a_tmp - 3.0) + 0.5 * K * ((s - 1.0) * (s - 1.0));
  e_a_tmp = K * (s - 1.0);
  for (jA = 0; jA < 6; jA++) {
    d = I3E[jA];
    d1 = a_tmp * (double)b[jA] - c_a_tmp * d;
    J1E[jA] = d1;
    d *= d_a_tmp;
    I3E[jA] = d;
    Se[jA] = C100 * d1 + e_a_tmp * d;
  }
  S[0] = Se[0];
  S[3] = Se[3];
  S[6] = Se[5];
  S[1] = Se[3];
  S[4] = Se[1];
  S[7] = Se[4];
  S[2] = Se[5];
  S[5] = Se[4];
  S[8] = Se[2];
  I3EE[0] = 0.0;
  I3EE[6] = 4.0 * C[8];
  I3EE[12] = 4.0 * C[4];
  I3EE[18] = 0.0;
  I3EE[24] = -4.0 * C[7];
  I3EE[30] = 0.0;
  I3EE[1] = 4.0 * C[8];
  I3EE[7] = 0.0;
  I3EE[13] = 4.0 * C[0];
  I3EE[19] = 0.0;
  I3EE[25] = 0.0;
  I3EE[31] = -4.0 * C[6];
  I3EE[2] = 4.0 * C[4];
  I3EE[8] = 4.0 * C[0];
  I3EE[14] = 0.0;
  I3EE[20] = -4.0 * C[3];
  I3EE[26] = 0.0;
  I3EE[32] = 0.0;
  I3EE[3] = 0.0;
  I3EE[9] = 0.0;
  I3EE[15] = -4.0 * C[3];
  I3EE[21] = -2.0 * C[8];
  I3EE[27] = 2.0 * C[6];
  I3EE[33] = 2.0 * C[7];
  I3EE[4] = -4.0 * C[7];
  I3EE[10] = 0.0;
  I3EE[16] = 0.0;
  I3EE[22] = 2.0 * C[6];
  I3EE[28] = -2.0 * C[0];
  I3EE[34] = 2.0 * C[3];
  I3EE[5] = 0.0;
  I3EE[11] = -4.0 * C[6];
  I3EE[17] = 0.0;
  I3EE[23] = 2.0 * C[7];
  I3EE[29] = 2.0 * C[3];
  I3EE[35] = -2.0 * C[4];
  /*  */
  smax = 0.88888888888888884 * I1 * b_a_tmp;
  /*  */
  s = -(0.66666666666666663 * a_tmp_tmp);
  /*  */
  for (i = 0; i < 6; i++) {
    for (jp1j = 0; jp1j < 6; jp1j++) {
      jA = jp1j + 6 * i;
      J1EE_tmp[jA] = I3E[jp1j] * I3E[i];
      D[jA] = J1E[jp1j] * I3E[i];
      b_I3E[jA] = I3E[jp1j] * J1E[i];
    }
  }
  for (i = 0; i < 36; i++) {
    d = J1EE_tmp[i];
    d1 = I3EE[i];
    D[i] =
        (C100 * ((s * (D[i] + b_I3E[i]) + smax * d) - c_a_tmp * d1) + K * d) +
        e_a_tmp * (-a_tmp_tmp * d + d_a_tmp * d1);
  }
}

/* End of code generation (PiollaStressNHFast.c) */
