/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgerc.c
 *
 * Code generation for function 'xgerc'
 *
 */

/* Include files */
#include "xgerc.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void xgerc(int m, int n, double alpha1, int ix0, const double y[3], double A[9],
           int ia0)
{
  double temp;
  int i;
  int ijA;
  int j;
  int jA;
  if (!(alpha1 == 0.0)) {
    jA = ia0;
    for (j = 0; j < n; j++) {
      temp = y[j];
      if (temp != 0.0) {
        temp *= alpha1;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }
      jA += 3;
    }
  }
}

/* End of code generation (xgerc.c) */
