/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * NonlinearStrainOperatorFast.c
 *
 * Code generation for function 'NonlinearStrainOperatorFast'
 *
 */

/* Include files */
#include "NonlinearStrainOperatorFast.h"
#include "NonlinearStrainOperatorFast_emxutil.h"
#include "NonlinearStrainOperatorFast_types.h"
#include <math.h>

/* Function Definitions */
void NonlinearStrainOperatorFast(const emxArray_real_T *N,
                                 const emxArray_real_T *dNdx, const double F[9],
                                 emxArray_real_T *Bn, emxArray_real_T *Bg,
                                 emxArray_real_T *NN, double *tau)
{
  emxArray_real_T *id1;
  emxArray_real_T *id2;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  double zz;
  int i;
  int i1;
  int i2;
  int loop_ub;
  zz = (double)dNdx->size[1] * (double)N->size[0];
  i = NN->size[0] * NN->size[1];
  NN->size[0] = dNdx->size[1];
  NN->size[1] = (int)zz;
  emxEnsureCapacity_real_T(NN, i);
  loop_ub = dNdx->size[1] * (int)zz;
  for (i = 0; i < loop_ub; i++) {
    NN->data[i] = 0.0;
  }
  emxInit_real_T(&id1, 2);
  if ((dNdx->size[1] == 0) || (zz < 1.0)) {
    id1->size[0] = 1;
    id1->size[1] = 0;
  } else {
    i = id1->size[0] * id1->size[1];
    id1->size[0] = 1;
    id1->size[1] = (int)floor((zz - 1.0) / (double)dNdx->size[1]) + 1;
    emxEnsureCapacity_real_T(id1, i);
    loop_ub = (int)floor((zz - 1.0) / (double)dNdx->size[1]);
    for (i = 0; i <= loop_ub; i++) {
      id1->data[i] = (double)dNdx->size[1] * (double)i + 1.0;
    }
  }
  emxInit_real_T(&id2, 2);
  if ((dNdx->size[1] == 0) || (zz < 2.0)) {
    id2->size[0] = 1;
    id2->size[1] = 0;
  } else {
    i = id2->size[0] * id2->size[1];
    id2->size[0] = 1;
    id2->size[1] = (int)floor((zz - 2.0) / (double)dNdx->size[1]) + 1;
    emxEnsureCapacity_real_T(id2, i);
    loop_ub = (int)floor((zz - 2.0) / (double)dNdx->size[1]);
    for (i = 0; i <= loop_ub; i++) {
      id2->data[i] = (double)dNdx->size[1] * (double)i + 2.0;
    }
  }
  emxInit_real_T(&r, 1);
  i = r->size[0];
  r->size[0] = id1->size[1];
  emxEnsureCapacity_real_T(r, i);
  loop_ub = id1->size[1];
  for (i = 0; i < loop_ub; i++) {
    r->data[i] = id1->data[i];
  }
  loop_ub = N->size[0];
  for (i = 0; i < loop_ub; i++) {
    NN->data[NN->size[0] * ((int)r->data[i] - 1)] = N->data[i];
  }
  emxInit_real_T(&r1, 1);
  i = r1->size[0];
  r1->size[0] = id2->size[1];
  emxEnsureCapacity_real_T(r1, i);
  loop_ub = id2->size[1];
  for (i = 0; i < loop_ub; i++) {
    r1->data[i] = id2->data[i];
  }
  emxFree_real_T(&id2);
  loop_ub = N->size[0];
  for (i = 0; i < loop_ub; i++) {
    NN->data[NN->size[0] * ((int)r1->data[i] - 1) + 1] = N->data[i];
  }
  if (dNdx->size[1] == 3) {
    if (3.0 > zz) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = 3;
    }
    loop_ub = N->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      NN->data[NN->size[0] * (i + i1 * i2) + 2] = N->data[i2];
    }
    loop_ub = dNdx->size[0];
    i = id1->size[0] * id1->size[1];
    id1->size[0] = 1;
    id1->size[1] = dNdx->size[0];
    emxEnsureCapacity_real_T(id1, i);
    for (i = 0; i < loop_ub; i++) {
      id1->data[i] = dNdx->data[i + dNdx->size[0] * 2];
    }
  } else {
    i = id1->size[0] * id1->size[1];
    id1->size[0] = 1;
    id1->size[1] = 1;
    emxEnsureCapacity_real_T(id1, i);
    id1->data[0] = 0.0;
  }
  i = Bn->size[0] * Bn->size[1];
  Bn->size[0] = (dNdx->size[1] - 1) * 3;
  Bn->size[1] = (int)zz;
  emxEnsureCapacity_real_T(Bn, i);
  loop_ub = (dNdx->size[1] - 1) * 3 * (int)zz;
  for (i = 0; i < loop_ub; i++) {
    Bn->data[i] = 0.0;
  }
  i = Bg->size[0] * Bg->size[1];
  Bg->size[0] = (int)(((double)dNdx->size[1] - 1.0) * 4.0 +
                      ((double)dNdx->size[1] - 2.0));
  Bg->size[1] = (int)zz;
  emxEnsureCapacity_real_T(Bg, i);
  loop_ub = (int)(((double)dNdx->size[1] - 1.0) * 4.0 +
                  ((double)dNdx->size[1] - 2.0)) *
            (int)zz;
  for (i = 0; i < loop_ub; i++) {
    Bg->data[i] = 0.0;
  }
  if (dNdx->size[1] == 2) {
    /*  2-dimensional */
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bn->data[Bn->size[0] * ((int)r->data[i] - 1)] = dNdx->data[i] * F[0];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bn->data[Bn->size[0] * ((int)r1->data[i] - 1)] = dNdx->data[i] * F[1];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bn->data[Bn->size[0] * ((int)r->data[i] - 1) + 1] =
          dNdx->data[i + dNdx->size[0]] * F[3];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bn->data[Bn->size[0] * ((int)r1->data[i] - 1) + 1] =
          dNdx->data[i + dNdx->size[0]] * F[4];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bn->data[Bn->size[0] * ((int)r->data[i] - 1) + 2] =
          dNdx->data[i] * F[3] + dNdx->data[i + dNdx->size[0]] * F[0];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bn->data[Bn->size[0] * ((int)r1->data[i] - 1) + 2] =
          dNdx->data[i] * F[4] + dNdx->data[i + dNdx->size[0]] * F[1];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bg->data[Bg->size[0] * ((int)r->data[i] - 1)] = dNdx->data[i];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bg->data[Bg->size[0] * ((int)r->data[i] - 1) + 1] =
          dNdx->data[i + dNdx->size[0]];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bg->data[Bg->size[0] * ((int)r1->data[i] - 1) + 2] = dNdx->data[i];
    }
    loop_ub = dNdx->size[0];
    for (i = 0; i < loop_ub; i++) {
      Bg->data[Bg->size[0] * ((int)r1->data[i] - 1) + 3] =
          dNdx->data[i + dNdx->size[0]];
    }
  } else {
    /*  3-dimensional */
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bn->data[Bn->size[0] * (i * i1)] = dNdx->data[i1] * F[0];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2)] = dNdx->data[i2] * F[1];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2)] = dNdx->data[i2] * F[2];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bn->data[Bn->size[0] * (i * i1) + 1] =
          dNdx->data[i1 + dNdx->size[0]] * F[3];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 1] =
          dNdx->data[i2 + dNdx->size[0]] * F[4];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 1] =
          dNdx->data[i2 + dNdx->size[0]] * F[5];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = id1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bn->data[Bn->size[0] * (i * i1) + 2] = id1->data[i1] * F[6];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = id1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 2] = id1->data[i2] * F[7];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = id1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 2] = id1->data[i2] * F[8];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bn->data[Bn->size[0] * (i * i1) + 3] =
          dNdx->data[i1] * F[3] + dNdx->data[i1 + dNdx->size[0]] * F[0];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 3] =
          dNdx->data[i2] * F[4] + dNdx->data[i2 + dNdx->size[0]] * F[1];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 3] =
          dNdx->data[i2] * F[5] + dNdx->data[i2 + dNdx->size[0]] * F[2];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bn->data[Bn->size[0] * (i * i1) + 4] =
          dNdx->data[i1 + dNdx->size[0]] * F[6] + id1->data[i1] * F[3];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 4] =
          dNdx->data[i2 + dNdx->size[0]] * F[7] + id1->data[i2] * F[4];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 4] =
          dNdx->data[i2 + dNdx->size[0]] * F[8] + id1->data[i2] * F[5];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bn->data[Bn->size[0] * (i * i1) + 5] =
          dNdx->data[i1] * F[6] + id1->data[i1] * F[0];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 5] =
          dNdx->data[i2] * F[7] + id1->data[i2] * F[1];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bn->data[Bn->size[0] * (i + i1 * i2) + 5] =
          dNdx->data[i2] * F[8] + id1->data[i2] * F[2];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bg->data[Bg->size[0] * (i * i1)] = dNdx->data[i1];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bg->data[Bg->size[0] * (i * i1) + 1] = dNdx->data[i1 + dNdx->size[0]];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (1.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 1;
    } else {
      i = dNdx->size[1];
    }
    loop_ub = id1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Bg->data[Bg->size[0] * (i * i1) + 2] = id1->data[i1];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bg->data[Bg->size[0] * (i + i1 * i2) + 3] = dNdx->data[i2];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bg->data[Bg->size[0] * (i + i1 * i2) + 4] =
          dNdx->data[i2 + dNdx->size[0]];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (2.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 1;
      i1 = dNdx->size[1];
    }
    loop_ub = id1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bg->data[Bg->size[0] * (i + i1 * i2) + 5] = id1->data[i2];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bg->data[Bg->size[0] * (i + i1 * i2) + 6] = dNdx->data[i2];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = dNdx->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bg->data[Bg->size[0] * (i + i1 * i2) + 7] =
          dNdx->data[i2 + dNdx->size[0]];
    }
    if ((dNdx->size[1] == 0) ||
        ((dNdx->size[1] > 0) &&
         (3.0 > (double)dNdx->size[1] * (double)N->size[0]))) {
      i = 0;
      i1 = 1;
    } else {
      i = 2;
      i1 = dNdx->size[1];
    }
    loop_ub = id1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      Bg->data[Bg->size[0] * (i + i1 * i2) + 8] = id1->data[i2];
    }
  }
  emxFree_real_T(&r1);
  emxFree_real_T(&r);
  emxFree_real_T(&id1);
  *tau = 1.0;
}

/* End of code generation (NonlinearStrainOperatorFast.c) */
