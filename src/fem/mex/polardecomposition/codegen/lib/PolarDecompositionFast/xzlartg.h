/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlartg.h
 *
 * Code generation for function 'xzlartg'
 *
 */

#ifndef XZLARTG_H
#define XZLARTG_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);

void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
             creal_T *r);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xzlartg.h) */
