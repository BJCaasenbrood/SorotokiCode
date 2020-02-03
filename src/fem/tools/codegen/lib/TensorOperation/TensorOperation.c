/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * TensorOperation.c
 *
 * Code generation for function 'TensorOperation'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "TensorOperation.h"

/* Function Definitions */
void TensorOperation(const double A[9], const double B[9], boolean_T Arg, double
                     T[36])
{
  if (Arg) {
    T[0] = A[0] * B[0];
    T[6] = A[4] * B[0];
    T[12] = A[8] * B[0];
    T[18] = 2.0 * A[3] * B[0];
    T[24] = 2.0 * A[7] * B[0];
    T[30] = 2.0 * A[6] * B[0];
    T[1] = A[0] * B[4];
    T[7] = A[4] * B[4];
    T[13] = A[8] * B[4];
    T[19] = 2.0 * A[3] * B[4];
    T[25] = 2.0 * A[7] * B[4];
    T[31] = 2.0 * A[6] * B[4];
    T[2] = A[0] * B[8];
    T[8] = A[4] * B[8];
    T[14] = A[8] * B[8];
    T[20] = 2.0 * A[3] * B[8];
    T[26] = 2.0 * A[7] * B[8];
    T[32] = 2.0 * A[6] * B[8];
    T[3] = 2.0 * A[0] * B[3];
    T[9] = 2.0 * A[4] * B[3];
    T[15] = 2.0 * A[8] * B[3];
    T[21] = 4.0 * A[3] * B[3];
    T[27] = 4.0 * A[7] * B[3];
    T[33] = 4.0 * A[6] * B[3];
    T[4] = 2.0 * A[0] * B[7];
    T[10] = 2.0 * A[4] * B[7];
    T[16] = 2.0 * A[8] * B[7];
    T[22] = 4.0 * A[3] * B[7];
    T[28] = 4.0 * A[7] * B[7];
    T[34] = 4.0 * A[6] * B[7];
    T[5] = 2.0 * A[0] * B[6];
    T[11] = 2.0 * A[4] * B[6];
    T[17] = 2.0 * A[8] * B[6];
    T[23] = 4.0 * A[3] * B[6];
    T[29] = 4.0 * A[7] * B[6];
    T[35] = 4.0 * A[6] * B[6];
  } else {
    T[0] = A[0] * B[0];
    T[6] = A[1] * B[1];
    T[12] = A[2] * B[2];
    T[18] = 2.0 * A[0] * B[1];
    T[24] = 2.0 * A[1] * B[2];
    T[30] = 2.0 * A[0] * B[2];
    T[1] = A[3] * B[3];
    T[7] = A[4] * B[4];
    T[13] = A[5] * B[5];
    T[19] = 2.0 * A[3] * B[4];
    T[25] = 2.0 * A[4] * B[5];
    T[31] = 2.0 * A[3] * B[5];
    T[2] = A[6] * B[6];
    T[8] = A[7] * B[7];
    T[14] = A[8] * B[8];
    T[20] = 2.0 * A[6] * B[7];
    T[26] = 2.0 * A[7] * B[8];
    T[32] = 2.0 * A[6] * B[8];
    T[3] = A[0] * B[3] + A[3] * B[0];
    T[9] = A[1] * B[4] + A[4] * B[1];
    T[15] = A[2] * B[5] + A[5] * B[2];
    T[21] = 2.0 * A[0] * B[4] + 2.0 * A[3] * B[1];
    T[27] = 2.0 * A[1] * B[5] + 2.0 * A[4] * B[2];
    T[33] = 2.0 * A[0] * B[5] + 2.0 * A[3] * B[2];
    T[4] = A[3] * B[6] + A[6] * B[3];
    T[10] = A[4] * B[7] + A[7] * B[4];
    T[16] = A[5] * B[8] + A[8] * B[5];
    T[22] = 2.0 * A[3] * B[7] + 2.0 * A[6] * B[4];
    T[28] = 2.0 * A[4] * B[8] + 2.0 * A[7] * B[5];
    T[34] = 2.0 * A[3] * B[8] + 2.0 * A[6] * B[5];
    T[5] = A[0] * B[6] + A[6] * B[0];
    T[11] = A[1] * B[7] + A[7] * B[1];
    T[17] = A[2] * B[8] + A[8] * B[2];
    T[23] = 2.0 * A[0] * B[7] + 2.0 * A[6] * B[1];
    T[29] = 2.0 * A[1] * B[8] + 2.0 * A[7] * B[2];
    T[35] = 2.0 * A[0] * B[8] + 2.0 * A[6] * B[2];
  }
}

/* End of code generation (TensorOperation.c) */
