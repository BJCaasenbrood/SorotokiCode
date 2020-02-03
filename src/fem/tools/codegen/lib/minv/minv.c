/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minv.c
 *
 * Code generation for function 'minv'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "minv.h"

/* Function Definitions */
void minv(const double A[9], double X[9])
{
  int i0;
  int k;
  signed char I[9];
  int j;
  double augmat[18];
  double b_augmat;
  int b_j;
  double c_augmat[6];
  int g;

  /*  A_inverse as the name suggests is the inverse of A */
  /*  A is a square matrix which is invertible */
  /*  if A is not a square matrix or A is square matrix but not invertible,then the output is not equal to inverse of A */
  /*  the order of A is axa */
  for (i0 = 0; i0 < 9; i0++) {
    I[i0] = 0;
  }

  for (k = 0; k < 3; k++) {
    I[k + 3 * k] = 1;
    for (i0 = 0; i0 < 3; i0++) {
      augmat[i0 + 3 * k] = A[i0 + 3 * k];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    for (j = 0; j < 3; j++) {
      augmat[j + 3 * (i0 + 3)] = I[j + 3 * i0];
    }
  }

  /*  AUGMENTED MATRIX */
  /*  GAUSSIAN ELMINATION METHOD: */
  /*  when A is invertible, [A I] is row equivalent to [I inv(A)] */
  /*  in other words, the row operations that convert A to I also converts I to inv(A) */
  /*  I is identity matrix of order axa, inv(A) is the inverse of A of order axa */
  /*  Converting A to its Echelon form */
  for (k = 0; k < 2; k++) {
    b_augmat = augmat[k + 3 * k];
    for (i0 = 0; i0 < 6; i0++) {
      augmat[k + 3 * i0] /= b_augmat;
    }

    /*  normalization,so that pivot values will be equal to 1 */
    for (j = 0; j <= 1 - k; j++) {
      b_j = (k + j) + 1;
      for (i0 = 0; i0 < 6; i0++) {
        c_augmat[i0] = augmat[b_j + 3 * i0] - augmat[k + 3 * i0] * augmat[b_j +
          3 * k];
      }

      for (i0 = 0; i0 < 6; i0++) {
        augmat[b_j + 3 * i0] = c_augmat[i0];
      }

      /*  making the elements below the pivots as zero */
    }
  }

  b_augmat = augmat[8];
  for (i0 = 0; i0 < 6; i0++) {
    augmat[2 + 3 * i0] /= b_augmat;
  }

  /*  noralization of the last row of A */
  /*  Converting A from its Echelon form to Row Reduced Echelon form */
  for (k = 0; k < 2; k++) {
    i0 = (int)((1.0 + (-1.0 - ((2.0 + (double)k) - 1.0))) / -1.0);
    for (b_j = 0; b_j < i0; b_j++) {
      g = k - b_j;
      for (j = 0; j < 6; j++) {
        c_augmat[j] = augmat[g + 3 * j] - augmat[(k + 3 * j) + 1] * augmat[g + 3
          * (k + 1)];
      }

      for (j = 0; j < 6; j++) {
        augmat[g + 3 * j] = c_augmat[j];
      }

      /*  makes the elements above pivots to be row */
    }
  }

  /* We end up with A converted to I and I will be converted to inv(A) */
  for (i0 = 0; i0 < 3; i0++) {
    for (j = 0; j < 3; j++) {
      X[j + 3 * i0] = augmat[j + 3 * (3 + i0)];
    }
  }

  /*  extracting inv(A) from augmented matrix [I inv(A)] */
}

/* End of code generation (minv.c) */
