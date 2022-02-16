/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * computeForwardKinematicsFast.c
 *
 * Code generation for function 'computeForwardKinematicsFast'
 *
 */

/* Include files */
#include "computeForwardKinematicsFast.h"
#include "computeForwardKinematicsFast_emxutil.h"
#include "computeForwardKinematicsFast_types.h"
#include <string.h>

/* Function Declarations */
static void ForwardODEX(const emxArray_real_T *x, const emxArray_real_T *Z1,
                        const emxArray_real_T *Theta, const double xia0[6],
                        const emxArray_real_T *Ba, emxArray_real_T *dZ1);

/* Function Definitions */
static void ForwardODEX(const emxArray_real_T *x, const emxArray_real_T *Z1,
                        const emxArray_real_T *Theta, const double xia0[6],
                        const emxArray_real_T *Ba, emxArray_real_T *dZ1)
{
  emxArray_real_T *BTh;
  emxArray_real_T *C;
  double A[36];
  double dv[9];
  double XI[6];
  double bkj;
  double d;
  double d1;
  int aoffset;
  int b_i;
  int boffset;
  int coffset;
  int i;
  int inner;
  int j;
  int k;
  int nc;
  emxInit_real_T(&BTh, 2);
  /* --------------------------------------------------------------------------
   */
  inner = Ba->size[1];
  nc = Theta->size[1];
  i = BTh->size[0] * BTh->size[1];
  BTh->size[0] = 6;
  BTh->size[1] = Theta->size[1];
  emxEnsureCapacity_real_T(BTh, i);
  for (j = 0; j < nc; j++) {
    coffset = j * 6;
    boffset = j * Theta->size[0];
    for (b_i = 0; b_i < 6; b_i++) {
      BTh->data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < inner; k++) {
      aoffset = k * 6;
      bkj = Theta->data[boffset + k];
      for (b_i = 0; b_i < 6; b_i++) {
        i = coffset + b_i;
        BTh->data[i] += Ba->data[aoffset + b_i] * bkj;
      }
    }
  }
  inner = BTh->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    XI[b_i] = 0.0;
  }
  for (k = 0; k < inner; k++) {
    aoffset = k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      XI[b_i] += BTh->data[aoffset + b_i] * x->data[k];
    }
  }
  for (i = 0; i < 6; i++) {
    XI[i] += xia0[i];
  }
  /*  build forward kin - position */
  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  memset(&A[0], 0, 36U * sizeof(double));
  for (i = 0; i < 3; i++) {
    bkj = Z1->data[6 * i];
    A[6 * i] = bkj;
    inner = 6 * (i + 3);
    A[inner + 3] = bkj;
    coffset = 6 * i + 1;
    bkj = Z1->data[coffset];
    A[coffset] = bkj;
    A[inner + 4] = bkj;
    coffset = 6 * i + 2;
    bkj = Z1->data[coffset];
    A[coffset] = bkj;
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
  for (i = 0; i < 3; i++) {
    bkj = dv[i];
    d = dv[i + 3];
    d1 = dv[i + 6];
    for (coffset = 0; coffset < 3; coffset++) {
      A[(i + 6 * coffset) + 3] =
          (bkj * Z1->data[6 * coffset] + d * Z1->data[6 * coffset + 1]) +
          d1 * Z1->data[6 * coffset + 2];
    }
  }
  i = dZ1->size[0] * dZ1->size[1];
  dZ1->size[0] = 6;
  dZ1->size[1] = (int)(((double)x->size[0] + 5.0) - 1.0);
  emxEnsureCapacity_real_T(dZ1, i);
  inner = 6 * (int)(((double)x->size[0] + 5.0) - 1.0);
  for (i = 0; i < inner; i++) {
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
    for (coffset = 0; coffset < 3; coffset++) {
      nc = i + 6 * coffset;
      dZ1->data[nc] = (Z1->data[i] * dv[3 * coffset] +
                       Z1->data[i + 6] * dv[3 * coffset + 1]) +
                      Z1->data[i + 12] * dv[3 * coffset + 2];
      dZ1->data[i + 18] += Z1->data[nc] * XI[coffset + 3];
    }
  }
  if (5U > x->size[0] + 4U) {
    i = 0;
  } else {
    i = 4;
  }
  emxInit_real_T(&C, 2);
  nc = BTh->size[1];
  coffset = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = BTh->size[1];
  emxEnsureCapacity_real_T(C, coffset);
  for (j = 0; j < nc; j++) {
    inner = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += A[k * 6 + b_i] * BTh->data[inner + k];
      }
      C->data[inner + b_i] = bkj;
    }
  }
  emxFree_real_T(&BTh);
  inner = C->size[1];
  for (coffset = 0; coffset < inner; coffset++) {
    for (nc = 0; nc < 6; nc++) {
      dZ1->data[nc + 6 * (i + coffset)] = C->data[nc + 6 * coffset];
    }
  }
  emxFree_real_T(&C);
}

void computeForwardKinematicsFast(const emxArray_real_T *x, double ds,
                                  const double p0[3], const double Phi0[9],
                                  const emxArray_real_T *xia0,
                                  const emxArray_real_T *Th,
                                  const emxArray_real_T *Ba,
                                  emxArray_real_T *gtmp, emxArray_real_T *Jtmp)
{
  emxArray_real_T *K1Z1;
  emxArray_real_T *K2Z1;
  emxArray_real_T *Z1;
  emxArray_real_T *b_Th;
  emxArray_real_T *b_Z1;
  double c_a[36];
  double Rt[9];
  double dv[9];
  double Ns;
  double a;
  double b_a;
  double d;
  double d1;
  int b_i;
  int b_loop_ub;
  int c_loop_ub;
  int d_loop_ub;
  int e_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int ii;
  int j;
  int k;
  int loop_ub;
  int n;
  emxInit_real_T(&Z1, 2);
  /*  % states */
  /*       % spatial steps */
  /*       % position zero */
  /*     % phi zeroclc */
  /*     % intrinsic strain vector */
  /*       % evaluated Theta matrix */
  /*  state to strain matrix         */
  /*  compute total length */
  Ns = (double)Th->size[2] / 2.0;
  i = Z1->size[0] * Z1->size[1];
  Z1->size[0] = 6;
  Z1->size[1] = (int)(((double)x->size[0] + 5.0) - 1.0);
  emxEnsureCapacity_real_T(Z1, i);
  loop_ub = 6 * (int)(((double)x->size[0] + 5.0) - 1.0);
  for (i = 0; i < loop_ub; i++) {
    Z1->data[i] = 0.0;
  }
  for (i = 0; i < 3; i++) {
    Z1->data[6 * i] = Phi0[3 * i];
    Z1->data[6 * i + 1] = Phi0[3 * i + 1];
    Z1->data[6 * i + 2] = Phi0[3 * i + 2];
    Z1->data[i + 18] = p0[i];
  }
  i = gtmp->size[0] * gtmp->size[1] * gtmp->size[2];
  gtmp->size[0] = 4;
  gtmp->size[1] = 4;
  i1 = (int)Ns;
  gtmp->size[2] = (int)Ns;
  emxEnsureCapacity_real_T(gtmp, i);
  i = Jtmp->size[0] * Jtmp->size[1] * Jtmp->size[2];
  Jtmp->size[0] = 6;
  Jtmp->size[1] = x->size[0];
  Jtmp->size[2] = (int)Ns;
  emxEnsureCapacity_real_T(Jtmp, i);
  loop_ub = 6 * x->size[0] * (int)Ns;
  for (i = 0; i < loop_ub; i++) {
    Jtmp->data[i] = 0.0;
  }
  if (0 <= (int)Ns - 1) {
    b_loop_ub = Th->size[0];
    i2 = Th->size[1];
    c_loop_ub = Th->size[1];
    a = 0.66666666666666663 * ds;
    d_loop_ub = Th->size[0];
    i3 = Th->size[1];
    e_loop_ub = Th->size[1];
    b_a = 0.25 * ds;
    if (5U > x->size[0] + 4U) {
      i4 = 0;
      i5 = -1;
    } else {
      i4 = 4;
      i5 = x->size[0] + 3;
    }
    dv[0] = 0.0;
    dv[4] = 0.0;
    dv[8] = 0.0;
    n = i5 - i4;
  }
  emxInit_real_T(&K1Z1, 2);
  emxInit_real_T(&K2Z1, 2);
  emxInit_real_T(&b_Th, 2);
  emxInit_real_T(&b_Z1, 2);
  for (ii = 0; ii < i1; ii++) {
    /*  first EL-diff eval */
    i = (ii + 1) << 1;
    j = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = b_loop_ub;
    b_Th->size[1] = i2;
    emxEnsureCapacity_real_T(b_Th, j);
    for (j = 0; j < c_loop_ub; j++) {
      for (b_i = 0; b_i < b_loop_ub; b_i++) {
        b_Th->data[b_i + b_Th->size[0] * j] =
            Th->data[(b_i + Th->size[0] * j) +
                     Th->size[0] * Th->size[1] * (i - 2)];
      }
    }
    ForwardODEX(x, Z1, b_Th, *(double(*)[6]) & xia0->data[6 * (i - 2)], Ba,
                K1Z1);
    /*  second EL-diff eval */
    j = b_Z1->size[0] * b_Z1->size[1];
    b_Z1->size[0] = 6;
    b_Z1->size[1] = Z1->size[1];
    emxEnsureCapacity_real_T(b_Z1, j);
    loop_ub = 6 * Z1->size[1];
    for (j = 0; j < loop_ub; j++) {
      b_Z1->data[j] = Z1->data[j] + a * K1Z1->data[j];
    }
    j = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = d_loop_ub;
    b_Th->size[1] = i3;
    emxEnsureCapacity_real_T(b_Th, j);
    for (j = 0; j < e_loop_ub; j++) {
      for (b_i = 0; b_i < d_loop_ub; b_i++) {
        b_Th->data[b_i + b_Th->size[0] * j] =
            Th->data[(b_i + Th->size[0] * j) +
                     Th->size[0] * Th->size[1] * (i - 1)];
      }
    }
    ForwardODEX(x, b_Z1, b_Th, *(double(*)[6]) & xia0->data[6 * (i - 1)], Ba,
                K2Z1);
    /*  update integrands */
    loop_ub = 6 * Z1->size[1];
    i = Z1->size[0] * Z1->size[1];
    Z1->size[0] = 6;
    emxEnsureCapacity_real_T(Z1, i);
    for (i = 0; i < loop_ub; i++) {
      Z1->data[i] += b_a * (K1Z1->data[i] + 3.0 * K2Z1->data[i]);
    }
    /*  recover the kinematics entities */
    /* --------------------------------------------------------------------------
     */
    for (i = 0; i < 4; i++) {
      j = 4 * i + 16 * ii;
      gtmp->data[j] = 0.0;
      gtmp->data[j + 1] = 0.0;
      gtmp->data[j + 2] = 0.0;
      gtmp->data[j + 3] = 0.0;
    }
    gtmp->data[16 * ii + 15] = 1.0;
    /* --------------------------------------------------------------------------
     */
    for (i = 0; i < 3; i++) {
      j = 4 * i + 16 * ii;
      gtmp->data[j] = Z1->data[6 * i];
      Rt[3 * i] = Z1->data[i];
      gtmp->data[j + 1] = Z1->data[6 * i + 1];
      Rt[3 * i + 1] = Z1->data[i + 6];
      gtmp->data[j + 2] = Z1->data[6 * i + 2];
      Rt[3 * i + 2] = Z1->data[i + 12];
      gtmp->data[(i + 16 * ii) + 12] = Z1->data[i + 18];
    }
    /* --------------------------------------------------------------------------
     */
    memset(&c_a[0], 0, 36U * sizeof(double));
    for (i = 0; i < 3; i++) {
      Ns = Rt[3 * i];
      c_a[6 * i] = Ns;
      loop_ub = 6 * (i + 3);
      c_a[loop_ub + 3] = Ns;
      Ns = Rt[3 * i + 1];
      c_a[6 * i + 1] = Ns;
      c_a[loop_ub + 4] = Ns;
      Ns = Rt[3 * i + 2];
      c_a[6 * i + 2] = Ns;
      c_a[loop_ub + 5] = Ns;
    }
    dv[1] = -Z1->data[20];
    dv[2] = Z1->data[19];
    dv[3] = Z1->data[20];
    dv[5] = -Z1->data[18];
    dv[6] = -Z1->data[19];
    dv[7] = Z1->data[18];
    for (i = 0; i < 3; i++) {
      Ns = Rt[i];
      d = Rt[i + 3];
      d1 = Rt[i + 6];
      for (j = 0; j < 3; j++) {
        c_a[(i + 6 * j) + 3] =
            (Ns * dv[3 * j] + d * dv[3 * j + 1]) + d1 * dv[3 * j + 2];
      }
    }
    i = K1Z1->size[0] * K1Z1->size[1];
    K1Z1->size[0] = 6;
    K1Z1->size[1] = (i5 - i4) + 1;
    emxEnsureCapacity_real_T(K1Z1, i);
    for (j = 0; j <= n; j++) {
      loop_ub = j * 6;
      for (b_i = 0; b_i < 6; b_i++) {
        Ns = 0.0;
        for (k = 0; k < 6; k++) {
          i = loop_ub + k;
          Ns += c_a[k * 6 + b_i] * Z1->data[i % 6 + 6 * (i4 + i / 6)];
        }
        K1Z1->data[loop_ub + b_i] = Ns;
      }
    }
    loop_ub = K1Z1->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (j = 0; j < 6; j++) {
        b_i = j + 6 * i;
        Jtmp->data[b_i + 6 * Jtmp->size[1] * ii] = K1Z1->data[b_i];
      }
    }
  }
  emxFree_real_T(&b_Z1);
  emxFree_real_T(&b_Th);
  emxFree_real_T(&K2Z1);
  emxFree_real_T(&K1Z1);
  emxFree_real_T(&Z1);
}

/* End of code generation (computeForwardKinematicsFast.c) */
