/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * computeLagrangianFast.c
 *
 * Code generation for function 'computeLagrangianFast'
 *
 */

/* Include files */
#include "computeLagrangianFast.h"
#include "computeLagrangianFast_emxutil.h"
#include "computeLagrangianFast_types.h"
#include "mtimes.h"
#include <string.h>

/* Function Declarations */
static void LagrangianODEX(const emxArray_real_T *x, const emxArray_real_T *dx,
                           const emxArray_real_T *Z1,
                           const emxArray_real_T *Theta, const double xia0[6],
                           const emxArray_real_T *Ba, const double Mtt[36],
                           const double Ktt[36], emxArray_real_T *dZ1,
                           emxArray_real_T *dZ2);

/* Function Definitions */
static void LagrangianODEX(const emxArray_real_T *x, const emxArray_real_T *dx,
                           const emxArray_real_T *Z1,
                           const emxArray_real_T *Theta, const double xia0[6],
                           const emxArray_real_T *Ba, const double Mtt[36],
                           const double Ktt[36], emxArray_real_T *dZ1,
                           emxArray_real_T *dZ2)
{
  static const short b[6] = {0, 0, 0, 0, 0, 9810};
  emxArray_real_T *C;
  emxArray_real_T *Jg;
  emxArray_real_T *b_C;
  emxArray_real_T *c_C;
  emxArray_real_T *dG;
  emxArray_real_T *d_C;
  emxArray_real_T *r;
  emxArray_real_T *y_tmp;
  double A[36];
  double Ai[36];
  double adV[36];
  double b_Ai[36];
  double Rt[9];
  double dv[9];
  double V[6];
  double XI[6];
  double y[6];
  double bkj;
  double d;
  double d1;
  int aoffset;
  int b_i;
  int boffset;
  int coffset;
  int i;
  int i1;
  int i2;
  int i3;
  int inner;
  int j;
  int k;
  int n;
  int nc;
  unsigned int u;
  if (5U > x->size[0] + 4U) {
    i = 0;
    i1 = -1;
  } else {
    i = 4;
    i1 = x->size[0] + 3;
  }
  bkj = 2.0 * ((double)x->size[0] - 1.0) + 6.0;
  if (((double)x->size[0] + 6.0) - 1.0 > bkj) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = x->size[0] + 4;
    i3 = (int)bkj;
  }
  emxInit_real_T(&y_tmp, 2);
  /* Theta_ = ThetaEval;%Model.ShpFnc(s); */
  inner = Ba->size[1];
  nc = Theta->size[1];
  n = y_tmp->size[0] * y_tmp->size[1];
  y_tmp->size[0] = 6;
  y_tmp->size[1] = Theta->size[1];
  emxEnsureCapacity_real_T(y_tmp, n);
  for (j = 0; j < nc; j++) {
    coffset = j * 6;
    boffset = j * Theta->size[0];
    for (b_i = 0; b_i < 6; b_i++) {
      y_tmp->data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < inner; k++) {
      aoffset = k * 6;
      bkj = Theta->data[boffset + k];
      for (b_i = 0; b_i < 6; b_i++) {
        n = coffset + b_i;
        y_tmp->data[n] += Ba->data[aoffset + b_i] * bkj;
      }
    }
  }
  inner = y_tmp->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    XI[b_i] = 0.0;
  }
  for (k = 0; k < inner; k++) {
    aoffset = k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      XI[b_i] += y_tmp->data[aoffset + b_i] * x->data[k];
    }
  }
  for (n = 0; n < 6; n++) {
    XI[n] += xia0[n];
  }
  /*  build forward kin - position */
  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  memset(&A[0], 0, 36U * sizeof(double));
  for (n = 0; n < 3; n++) {
    bkj = Z1->data[6 * n];
    A[6 * n] = bkj;
    inner = 6 * (n + 3);
    A[inner + 3] = bkj;
    nc = 6 * n + 1;
    bkj = Z1->data[nc];
    A[nc] = bkj;
    A[inner + 4] = bkj;
    nc = 6 * n + 2;
    bkj = Z1->data[nc];
    A[nc] = bkj;
    A[inner + 5] = bkj;
  }
  dv[0] = 0.0;
  dv[3] = -Z1->data[20];
  dv[6] = Z1->data[19];
  dv[1] = Z1->data[20];
  dv[4] = 0.0;
  dv[7] = -Z1->data[18];
  dv[2] = -Z1->data[19];
  dv[5] = Z1->data[18];
  dv[8] = 0.0;
  /* --------------------------------------------------------------------------
   */
  for (n = 0; n < 3; n++) {
    bkj = dv[n];
    d = dv[n + 3];
    d1 = dv[n + 6];
    for (nc = 0; nc < 3; nc++) {
      A[(n + 6 * nc) + 3] =
          (bkj * Z1->data[6 * nc] + d * Z1->data[6 * nc + 1]) +
          d1 * Z1->data[6 * nc + 2];
    }
    Rt[3 * n] = Z1->data[n];
    Rt[3 * n + 1] = Z1->data[n + 6];
    Rt[3 * n + 2] = Z1->data[n + 12];
  }
  /* --------------------------------------------------------------------------
   */
  memset(&Ai[0], 0, 36U * sizeof(double));
  for (n = 0; n < 3; n++) {
    bkj = Rt[3 * n];
    Ai[6 * n] = bkj;
    inner = 6 * (n + 3);
    Ai[inner + 3] = bkj;
    bkj = Rt[3 * n + 1];
    Ai[6 * n + 1] = bkj;
    Ai[inner + 4] = bkj;
    bkj = Rt[3 * n + 2];
    Ai[6 * n + 2] = bkj;
    Ai[inner + 5] = bkj;
  }
  dv[0] = 0.0;
  dv[1] = -Z1->data[20];
  dv[2] = Z1->data[19];
  dv[3] = Z1->data[20];
  dv[4] = 0.0;
  dv[5] = -Z1->data[18];
  dv[6] = -Z1->data[19];
  dv[7] = Z1->data[18];
  dv[8] = 0.0;
  for (n = 0; n < 3; n++) {
    bkj = Rt[n];
    d = Rt[n + 3];
    d1 = Rt[n + 6];
    for (nc = 0; nc < 3; nc++) {
      Ai[(n + 6 * nc) + 3] =
          (bkj * dv[3 * nc] + d * dv[3 * nc + 1]) + d1 * dv[3 * nc + 2];
    }
  }
  emxInit_real_T(&Jg, 2);
  /*  build jacobian */
  inner = i1 - i;
  i1 = Jg->size[0] * Jg->size[1];
  Jg->size[0] = 6;
  Jg->size[1] = inner + 1;
  emxEnsureCapacity_real_T(Jg, i1);
  for (j = 0; j <= inner; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        i1 = nc + k;
        bkj += Ai[k * 6 + b_i] * Z1->data[i1 % 6 + 6 * (i + i1 / 6)];
      }
      Jg->data[nc + b_i] = bkj;
    }
  }
  inner = Jg->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    V[b_i] = 0.0;
  }
  for (k = 0; k < inner; k++) {
    aoffset = k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      V[b_i] += Jg->data[aoffset + b_i] * dx->data[k];
    }
  }
  /* --------------------------------------------------------------------------
   */
  memset(&adV[0], 0, 36U * sizeof(double));
  /* --------------------------------------------------------------------------
   */
  Rt[0] = 0.0;
  Rt[3] = -V[2];
  Rt[6] = V[1];
  Rt[1] = V[2];
  Rt[4] = 0.0;
  Rt[7] = -V[0];
  Rt[2] = -V[1];
  Rt[5] = V[0];
  Rt[8] = 0.0;
  /* --------------------------------------------------------------------------
   */
  for (i = 0; i < 3; i++) {
    bkj = Rt[3 * i];
    adV[6 * i] = bkj;
    inner = 6 * (i + 3);
    adV[inner + 3] = bkj;
    bkj = Rt[3 * i + 1];
    adV[6 * i + 1] = bkj;
    adV[inner + 4] = bkj;
    bkj = Rt[3 * i + 2];
    adV[6 * i + 2] = bkj;
    adV[inner + 5] = bkj;
  }
  adV[3] = 0.0;
  adV[9] = -V[5];
  adV[15] = V[4];
  adV[4] = V[5];
  adV[10] = 0.0;
  adV[16] = -V[3];
  adV[5] = -V[4];
  adV[11] = V[3];
  adV[17] = 0.0;
  /*  compute inertia, coriolis, gravity */
  for (i = 0; i < 6; i++) {
    y[i] = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      bkj = 0.0;
      for (n = 0; n < 6; n++) {
        bkj += Ai[i + 6 * n] * Mtt[n + 6 * i1];
      }
      y[i] += bkj * (double)b[i1];
    }
  }
  emxInit_real_T(&dG, 1);
  inner = Jg->size[1];
  i = dG->size[0];
  dG->size[0] = Jg->size[1];
  emxEnsureCapacity_real_T(dG, i);
  for (b_i = 0; b_i < inner; b_i++) {
    aoffset = b_i * 6;
    bkj = 0.0;
    for (k = 0; k < 6; k++) {
      bkj += Jg->data[aoffset + k] * y[k];
    }
    dG->data[b_i] = bkj;
  }
  /*  compute (nonlinear stiffness) */
  /*  compute grav. potential energy */
  i = dZ1->size[0] * dZ1->size[1];
  dZ1->size[0] = 6;
  dZ1->size[1] = (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  emxEnsureCapacity_real_T(dZ1, i);
  nc = 6 * (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  for (i = 0; i < nc; i++) {
    dZ1->data[i] = 0.0;
  }
  dv[0] = 0.0;
  dv[3] = -XI[2];
  dv[6] = XI[1];
  dv[1] = XI[2];
  dv[4] = 0.0;
  dv[7] = -XI[0];
  dv[2] = -XI[1];
  dv[5] = XI[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    dZ1->data[i + 18] = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      n = i + 6 * i1;
      dZ1->data[n] =
          (Z1->data[i] * dv[3 * i1] + Z1->data[i + 6] * dv[3 * i1 + 1]) +
          Z1->data[i + 12] * dv[3 * i1 + 2];
      dZ1->data[i + 18] += Z1->data[n] * XI[i1 + 3];
    }
  }
  if (5U > x->size[0] + 4U) {
    i = 0;
  } else {
    i = 4;
  }
  emxInit_real_T(&C, 2);
  n = y_tmp->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = y_tmp->size[1];
  emxEnsureCapacity_real_T(C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += A[k * 6 + b_i] * y_tmp->data[nc + k];
      }
      C->data[nc + b_i] = bkj;
    }
  }
  nc = C->size[1];
  for (i1 = 0; i1 < nc; i1++) {
    for (n = 0; n < 6; n++) {
      dZ1->data[n + 6 * (i + i1)] = C->data[n + 6 * i1];
    }
  }
  if (((double)x->size[0] + 6.0) - 1.0 >
      2.0 * ((double)x->size[0] - 1.0) + 6.0) {
    i = -4;
  } else {
    i = x->size[0];
  }
  for (i1 = 0; i1 < 6; i1++) {
    for (n = 0; n < 6; n++) {
      bkj = 0.0;
      for (nc = 0; nc < 6; nc++) {
        bkj += A[i1 + 6 * nc] * adV[nc + 6 * n];
      }
      b_Ai[i1 + 6 * n] = bkj;
    }
  }
  n = y_tmp->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = y_tmp->size[1];
  emxEnsureCapacity_real_T(C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += b_Ai[k * 6 + b_i] * y_tmp->data[nc + k];
      }
      C->data[nc + b_i] = bkj;
    }
  }
  nc = C->size[1];
  for (i1 = 0; i1 < nc; i1++) {
    for (n = 0; n < 6; n++) {
      dZ1->data[n + 6 * ((i + i1) + 4)] = C->data[n + 6 * i1];
    }
  }
  dZ1->data[22] =
      (Mtt[21] * Z1->data[18] * 0.0 + Mtt[21] * Z1->data[19] * 0.0) +
      Mtt[21] * Z1->data[20] * 9810.0;
  bkj = 0.0;
  for (i = 0; i < 6; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      d += 0.5 * V[i1] * Mtt[i1 + 6 * i];
    }
    bkj += d * V[i];
  }
  dZ1->data[23] = bkj;
  i = dZ2->size[0] * dZ2->size[1];
  dZ2->size[0] = x->size[0];
  dZ2->size[1] = (int)(3.0 * (double)x->size[0] + 1.0);
  emxEnsureCapacity_real_T(dZ2, i);
  nc = x->size[0] * (int)(3.0 * (double)x->size[0] + 1.0);
  for (i = 0; i < nc; i++) {
    dZ2->data[i] = 0.0;
  }
  emxInit_real_T(&b_C, 2);
  emxInit_real_T(&r, 2);
  mtimes(Jg, Mtt, r);
  b_mtimes(r, Jg, b_C);
  nc = b_C->size[1];
  for (i = 0; i < nc; i++) {
    inner = b_C->size[0];
    for (i1 = 0; i1 < inner; i1++) {
      dZ2->data[i1 + dZ2->size[0] * i] = b_C->data[i1 + b_C->size[0] * i];
    }
  }
  if (x->size[0] + 1U > ((unsigned int)x->size[0] << 1)) {
    i = 0;
  } else {
    i = x->size[0];
  }
  for (i1 = 0; i1 < 6; i1++) {
    for (n = 0; n < 6; n++) {
      bkj = 0.0;
      d = 0.0;
      for (nc = 0; nc < 6; nc++) {
        inner = nc + 6 * n;
        bkj += Mtt[i1 + 6 * nc] * adV[inner];
        d += adV[nc + 6 * i1] * Mtt[inner];
      }
      inner = i1 + 6 * n;
      b_Ai[inner] = d;
      A[inner] = bkj;
    }
  }
  for (i1 = 0; i1 < 36; i1++) {
    A[i1] -= b_Ai[i1];
  }
  n = Jg->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = Jg->size[1];
  emxEnsureCapacity_real_T(C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += A[k * 6 + b_i] * Jg->data[nc + k];
      }
      C->data[nc + b_i] = bkj;
    }
  }
  emxInit_real_T(&c_C, 2);
  inner = i3 - i2;
  n = inner - 1;
  i1 = c_C->size[0] * c_C->size[1];
  c_C->size[0] = 6;
  c_C->size[1] = inner;
  emxEnsureCapacity_real_T(c_C, i1);
  for (j = 0; j <= n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        i1 = nc + k;
        bkj += Ai[k * 6 + b_i] * Z1->data[i1 % 6 + 6 * (i2 + i1 / 6)];
      }
      c_C->data[nc + b_i] = bkj;
    }
  }
  emxInit_real_T(&d_C, 2);
  n = c_C->size[1];
  i1 = d_C->size[0] * d_C->size[1];
  d_C->size[0] = 6;
  d_C->size[1] = c_C->size[1];
  emxEnsureCapacity_real_T(d_C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += Mtt[k * 6 + b_i] * c_C->data[nc + k];
      }
      d_C->data[nc + b_i] = bkj;
    }
  }
  emxFree_real_T(&c_C);
  nc = 6 * C->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  emxEnsureCapacity_real_T(C, i1);
  for (i1 = 0; i1 < nc; i1++) {
    C->data[i1] += d_C->data[i1];
  }
  emxFree_real_T(&d_C);
  inner = Jg->size[1];
  n = C->size[1];
  i1 = b_C->size[0] * b_C->size[1];
  b_C->size[0] = Jg->size[1];
  b_C->size[1] = C->size[1];
  emxEnsureCapacity_real_T(b_C, i1);
  for (j = 0; j < n; j++) {
    coffset = j * inner;
    boffset = j * 6;
    for (b_i = 0; b_i < inner; b_i++) {
      aoffset = b_i * 6;
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += Jg->data[aoffset + k] * C->data[boffset + k];
      }
      b_C->data[coffset + b_i] = bkj;
    }
  }
  emxFree_real_T(&C);
  emxFree_real_T(&Jg);
  nc = b_C->size[1];
  for (i1 = 0; i1 < nc; i1++) {
    inner = b_C->size[0];
    for (i2 = 0; i2 < inner; i2++) {
      dZ2->data[i2 + dZ2->size[0] * (i + i1)] =
          b_C->data[i2 + b_C->size[0] * i1];
    }
  }
  u = ((unsigned int)x->size[0] << 1) + 1U;
  if (u > 3.0 * (double)x->size[0]) {
    i = 1;
  } else {
    i = (int)u;
  }
  mtimes(y_tmp, Ktt, r);
  b_mtimes(r, y_tmp, b_C);
  nc = b_C->size[1];
  emxFree_real_T(&r);
  emxFree_real_T(&y_tmp);
  for (i1 = 0; i1 < nc; i1++) {
    inner = b_C->size[0];
    for (i2 = 0; i2 < inner; i2++) {
      dZ2->data[i2 + dZ2->size[0] * ((i + i1) - 1)] =
          b_C->data[i2 + b_C->size[0] * i1];
    }
  }
  emxFree_real_T(&b_C);
  i = (int)(3.0 * (double)x->size[0] + 1.0) - 1;
  nc = dG->size[0];
  for (i1 = 0; i1 < nc; i1++) {
    dZ2->data[i1 + dZ2->size[0] * i] = dG->data[i1];
  }
  emxFree_real_T(&dG);
}

void computeLagrangianFast(const emxArray_real_T *x, const emxArray_real_T *dx,
                           double ds, const double p0[3], const double Phi0[9],
                           const emxArray_real_T *xia0,
                           const emxArray_real_T *Th, const emxArray_real_T *Ba,
                           const double Ktt[36], const double Mtt[36],
                           double Zeta, emxArray_real_T *M, emxArray_real_T *C,
                           emxArray_real_T *K, emxArray_real_T *R,
                           emxArray_real_T *G, double p[3], double Phi[9],
                           emxArray_real_T *J, double *Vg, double *Kin)
{
  emxArray_real_T *K1Z1;
  emxArray_real_T *K1Z2;
  emxArray_real_T *K2Z1;
  emxArray_real_T *K2Z2;
  emxArray_real_T *Z1;
  emxArray_real_T *Z2;
  emxArray_real_T *b_Th;
  emxArray_real_T *b_Z1;
  double a[36];
  double Rt[9];
  double dv[9];
  double d;
  double d1;
  double s;
  int b_i;
  int b_loop_ub;
  int i;
  int i1;
  int ii;
  int k;
  int loop_ub;
  unsigned int u;
  emxInit_real_T(&Z1, 2);
  /*  % states */
  /*       % spatial steps */
  /*       % position zero */
  /*     % phi zero */
  /*     % intrinsic strain vector */
  /*       % evaluated Theta matrix */
  /*       % state to strain matrix */
  /*      % geometric stiffness */
  /*      % geometric inertia */
  /*  compute total length */
  i = Z1->size[0] * Z1->size[1];
  Z1->size[0] = 6;
  Z1->size[1] = (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  emxEnsureCapacity_real_T(Z1, i);
  loop_ub = 6 * (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  for (i = 0; i < loop_ub; i++) {
    Z1->data[i] = 0.0;
  }
  emxInit_real_T(&Z2, 2);
  i = Z2->size[0] * Z2->size[1];
  Z2->size[0] = x->size[0];
  Z2->size[1] = (int)(3.0 * (double)x->size[0] + 1.0);
  emxEnsureCapacity_real_T(Z2, i);
  loop_ub = x->size[0] * (int)(3.0 * (double)x->size[0] + 1.0);
  for (i = 0; i < loop_ub; i++) {
    Z2->data[i] = 0.0;
  }
  for (i = 0; i < 3; i++) {
    Z1->data[6 * i] = Phi0[3 * i];
    Z1->data[6 * i + 1] = Phi0[3 * i + 1];
    Z1->data[6 * i + 2] = Phi0[3 * i + 2];
    Z1->data[i + 18] = p0[i];
  }
  /* NLStiff = false;  */
  i = (int)((double)Th->size[2] / 2.0);
  emxInit_real_T(&K1Z1, 2);
  emxInit_real_T(&K1Z2, 2);
  emxInit_real_T(&K2Z1, 2);
  emxInit_real_T(&K2Z2, 2);
  emxInit_real_T(&b_Th, 2);
  emxInit_real_T(&b_Z1, 2);
  for (ii = 0; ii < i; ii++) {
    /*  first EL-diff eval */
    i1 = (int)((unsigned int)(ii + 1) << 1);
    loop_ub = Th->size[0];
    b_loop_ub = Th->size[1];
    k = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = Th->size[0];
    b_Th->size[1] = Th->size[1];
    emxEnsureCapacity_real_T(b_Th, k);
    for (k = 0; k < b_loop_ub; k++) {
      for (b_i = 0; b_i < loop_ub; b_i++) {
        b_Th->data[b_i + b_Th->size[0] * k] =
            Th->data[(b_i + Th->size[0] * k) +
                     Th->size[0] * Th->size[1] * (i1 - 2)];
      }
    }
    LagrangianODEX(x, dx, Z1, b_Th,
                   *(double(*)[6]) &
                       xia0->data[6 * ((int)((unsigned int)(ii + 1) << 1) - 2)],
                   Ba, Mtt, Ktt, K1Z1, K1Z2);
    /*  second EL-diff eval */
    s = 0.66666666666666663 * ds;
    k = b_Z1->size[0] * b_Z1->size[1];
    b_Z1->size[0] = 6;
    b_Z1->size[1] = Z1->size[1];
    emxEnsureCapacity_real_T(b_Z1, k);
    loop_ub = 6 * Z1->size[1];
    for (k = 0; k < loop_ub; k++) {
      b_Z1->data[k] = Z1->data[k] + s * K1Z1->data[k];
    }
    loop_ub = Th->size[0];
    b_loop_ub = Th->size[1];
    k = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = Th->size[0];
    b_Th->size[1] = Th->size[1];
    emxEnsureCapacity_real_T(b_Th, k);
    for (k = 0; k < b_loop_ub; k++) {
      for (b_i = 0; b_i < loop_ub; b_i++) {
        b_Th->data[b_i + b_Th->size[0] * k] =
            Th->data[(b_i + Th->size[0] * k) +
                     Th->size[0] * Th->size[1] * (i1 - 1)];
      }
    }
    LagrangianODEX(x, dx, b_Z1, b_Th,
                   *(double(*)[6]) &
                       xia0->data[6 * ((int)((unsigned int)(ii + 1) << 1) - 1)],
                   Ba, Mtt, Ktt, K2Z1, K2Z2);
    /*  update integrands */
    s = 0.25 * ds;
    loop_ub = 6 * Z1->size[1];
    i1 = Z1->size[0] * Z1->size[1];
    Z1->size[0] = 6;
    emxEnsureCapacity_real_T(Z1, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      Z1->data[i1] += s * (K1Z1->data[i1] + 3.0 * K2Z1->data[i1]);
    }
    loop_ub = Z2->size[0] * Z2->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Z2->data[i1] += s * (K1Z2->data[i1] + 3.0 * K2Z2->data[i1]);
    }
  }
  emxFree_real_T(&b_Z1);
  emxFree_real_T(&b_Th);
  emxFree_real_T(&K2Z2);
  emxFree_real_T(&K2Z1);
  emxFree_real_T(&K1Z2);
  emxFree_real_T(&K1Z1);
  /*  recover the kinematics entities */
  for (i = 0; i < 3; i++) {
    p[i] = Z1->data[i + 18];
    Phi[3 * i] = Z1->data[6 * i];
    Phi[3 * i + 1] = Z1->data[6 * i + 1];
    Phi[3 * i + 2] = Z1->data[6 * i + 2];
  }
  if (5U > x->size[0] + 4U) {
    i = 0;
    i1 = -1;
  } else {
    i = 4;
    i1 = x->size[0] + 3;
  }
  /* --------------------------------------------------------------------------
   */
  for (k = 0; k < 3; k++) {
    Rt[3 * k] = Z1->data[k];
    Rt[3 * k + 1] = Z1->data[k + 6];
    Rt[3 * k + 2] = Z1->data[k + 12];
  }
  /* --------------------------------------------------------------------------
   */
  memset(&a[0], 0, 36U * sizeof(double));
  for (k = 0; k < 3; k++) {
    s = Rt[3 * k];
    a[6 * k] = s;
    ii = 6 * (k + 3);
    a[ii + 3] = s;
    s = Rt[3 * k + 1];
    a[6 * k + 1] = s;
    a[ii + 4] = s;
    s = Rt[3 * k + 2];
    a[6 * k + 2] = s;
    a[ii + 5] = s;
  }
  dv[0] = 0.0;
  dv[1] = -Z1->data[20];
  dv[2] = Z1->data[19];
  dv[3] = Z1->data[20];
  dv[4] = 0.0;
  dv[5] = -Z1->data[18];
  dv[6] = -Z1->data[19];
  dv[7] = Z1->data[18];
  dv[8] = 0.0;
  for (k = 0; k < 3; k++) {
    s = Rt[k];
    d = Rt[k + 3];
    d1 = Rt[k + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      a[(k + 6 * b_i) + 3] =
          (s * dv[3 * b_i] + d * dv[3 * b_i + 1]) + d1 * dv[3 * b_i + 2];
    }
  }
  ii = i1 - i;
  i1 = J->size[0] * J->size[1];
  J->size[0] = 6;
  J->size[1] = ii + 1;
  emxEnsureCapacity_real_T(J, i1);
  for (b_loop_ub = 0; b_loop_ub <= ii; b_loop_ub++) {
    loop_ub = b_loop_ub * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      s = 0.0;
      for (k = 0; k < 6; k++) {
        i1 = loop_ub + k;
        s += a[k * 6 + b_i] * Z1->data[i1 % 6 + 6 * (i + i1 / 6)];
      }
      J->data[loop_ub + b_i] = s;
    }
  }
  /*  recover the dynamics entities */
  if (1 > x->size[0]) {
    loop_ub = 0;
    b_loop_ub = 0;
  } else {
    loop_ub = x->size[0];
    b_loop_ub = x->size[0];
  }
  i = M->size[0] * M->size[1];
  M->size[0] = loop_ub;
  M->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(M, i);
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      M->data[i1 + M->size[0] * i] = Z2->data[i1 + Z2->size[0] * i];
    }
  }
  if (1 > x->size[0]) {
    loop_ub = 0;
  } else {
    loop_ub = x->size[0];
  }
  u = (unsigned int)x->size[0] << 1;
  if (x->size[0] + 1U > u) {
    i = 0;
    i1 = 0;
  } else {
    i = x->size[0];
    i1 = (int)u;
  }
  k = C->size[0] * C->size[1];
  C->size[0] = loop_ub;
  b_loop_ub = i1 - i;
  C->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(C, k);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (k = 0; k < loop_ub; k++) {
      C->data[k + C->size[0] * i1] = Z2->data[k + Z2->size[0] * (i + i1)];
    }
  }
  if (1 > x->size[0]) {
    loop_ub = 0;
  } else {
    loop_ub = x->size[0];
  }
  u = ((unsigned int)x->size[0] << 1) + 1U;
  s = 3.0 * (double)x->size[0];
  if (u > s) {
    i = 0;
    i1 = 0;
  } else {
    i = (int)u - 1;
    i1 = (int)s;
  }
  k = K->size[0] * K->size[1];
  K->size[0] = loop_ub;
  b_loop_ub = i1 - i;
  K->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(K, k);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (k = 0; k < loop_ub; k++) {
      K->data[k + K->size[0] * i1] = Z2->data[k + Z2->size[0] * (i + i1)];
    }
  }
  if (1 > x->size[0]) {
    ii = 0;
  } else {
    ii = x->size[0];
  }
  i1 = (int)(3.0 * (double)x->size[0] + 1.0);
  k = G->size[0];
  G->size[0] = ii;
  emxEnsureCapacity_real_T(G, k);
  for (k = 0; k < ii; k++) {
    G->data[k] = Z2->data[k + Z2->size[0] * (i1 - 1)];
  }
  *Vg = Z1->data[22];
  *Kin = Z1->data[23];
  i1 = R->size[0] * R->size[1];
  R->size[0] = loop_ub;
  R->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(R, i1);
  emxFree_real_T(&Z1);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (k = 0; k < loop_ub; k++) {
      R->data[k + R->size[0] * i1] =
          Zeta * Z2->data[k + Z2->size[0] * (i + i1)];
    }
  }
  emxFree_real_T(&Z2);
}

/* End of code generation (computeLagrangianFast.c) */
