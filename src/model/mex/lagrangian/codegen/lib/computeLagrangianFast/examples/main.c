/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "computeLagrangianFast.h"
#include "computeLagrangianFast_emxAPI.h"
#include "computeLagrangianFast_terminate.h"
#include "computeLagrangianFast_types.h"

/* Function Declarations */
static void argInit_3x1_real_T(double result[3]);

static void argInit_3x3_real_T(double result[9]);

static emxArray_real_T *argInit_6x1xUnbounded_real_T(void);

static void argInit_6x6_real_T(double result[36]);

static emxArray_real_T *argInit_6xUnbounded_real_T(void);

static emxArray_real_T *argInit_Unboundedx1_real_T(void);

static double argInit_real_T(void);

static emxArray_real_T *c_argInit_UnboundedxUnboundedxU(void);

static void main_computeLagrangianFast(void);

/* Function Definitions */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

static void argInit_3x3_real_T(double result[9])
{
  int idx0;
  int idx1;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real_T();
    }
  }
}

static emxArray_real_T *argInit_6x1xUnbounded_real_T(void)
{
  emxArray_real_T *result;
  const int iv[3] = {6, 1, 2};
  int idx0;
  int idx1;
  int idx2;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result = emxCreateND_real_T(3, &iv[0]);
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 6; idx0++) {
    for (idx1 = 0; idx1 < 1; idx1++) {
      for (idx2 = 0; idx2 < result->size[2U]; idx2++) {
        /* Set the value of the array element.
Change this value to the value that the application requires. */
        result->data[idx0 + 6 * idx2] = argInit_real_T();
      }
    }
  }
  return result;
}

static void argInit_6x6_real_T(double result[36])
{
  int idx0;
  int idx1;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 6; idx0++) {
    for (idx1 = 0; idx1 < 6; idx1++) {
      /* Set the value of the array element.
Change this value to the value that the application requires. */
      result[idx0 + 6 * idx1] = argInit_real_T();
    }
  }
}

static emxArray_real_T *argInit_6xUnbounded_real_T(void)
{
  emxArray_real_T *result;
  int idx0;
  int idx1;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result = emxCreate_real_T(6, 2);
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 6; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
Change this value to the value that the application requires. */
      result->data[idx0 + 6 * idx1] = argInit_real_T();
    }
  }
  return result;
}

static emxArray_real_T *argInit_Unboundedx1_real_T(void)
{
  emxArray_real_T *result;
  const int i = 2;
  int idx0;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result = emxCreateND_real_T(1, &i);
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result->data[idx0] = argInit_real_T();
  }
  return result;
}

static double argInit_real_T(void)
{
  return 0.0;
}

static emxArray_real_T *c_argInit_UnboundedxUnboundedxU(void)
{
  emxArray_real_T *result;
  const int iv[3] = {2, 2, 2};
  int idx0;
  int idx1;
  int idx2;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result = emxCreateND_real_T(3, &iv[0]);
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      for (idx2 = 0; idx2 < result->size[2U]; idx2++) {
        /* Set the value of the array element.
Change this value to the value that the application requires. */
        result->data[(idx0 + result->size[0] * idx1) +
                     result->size[0] * result->size[1] * idx2] =
            argInit_real_T();
      }
    }
  }
  return result;
}

static void main_computeLagrangianFast(void)
{
  emxArray_real_T *Ba;
  emxArray_real_T *C;
  emxArray_real_T *G;
  emxArray_real_T *J;
  emxArray_real_T *K;
  emxArray_real_T *M;
  emxArray_real_T *R;
  emxArray_real_T *Th;
  emxArray_real_T *dx;
  emxArray_real_T *x;
  emxArray_real_T *xia0;
  double Ktt_tmp[36];
  double Phi[9];
  double dv[9];
  double p[3];
  double p0_tmp[3];
  double Kin;
  double Vg;
  double ds_tmp;
  emxInitArray_real_T(&M, 2);
  emxInitArray_real_T(&C, 2);
  emxInitArray_real_T(&K, 2);
  emxInitArray_real_T(&R, 2);
  emxInitArray_real_T(&G, 1);
  emxInitArray_real_T(&J, 2);
  /* Initialize function 'computeLagrangianFast' input arguments. */
  /* Initialize function input argument 'x'. */
  x = argInit_Unboundedx1_real_T();
  /* Initialize function input argument 'dx'. */
  dx = argInit_Unboundedx1_real_T();
  ds_tmp = argInit_real_T();
  /* Initialize function input argument 'p0'. */
  argInit_3x1_real_T(p0_tmp);
  /* Initialize function input argument 'Phi0'. */
  /* Initialize function input argument 'xia0'. */
  xia0 = argInit_6x1xUnbounded_real_T();
  /* Initialize function input argument 'Th'. */
  Th = c_argInit_UnboundedxUnboundedxU();
  /* Initialize function input argument 'Ba'. */
  Ba = argInit_6xUnbounded_real_T();
  /* Initialize function input argument 'Ktt'. */
  argInit_6x6_real_T(Ktt_tmp);
  /* Initialize function input argument 'Mtt'. */
  /* Initialize function input argument 'Gvec'. */
  /* Call the entry-point 'computeLagrangianFast'. */
  argInit_3x3_real_T(dv);
  computeLagrangianFast(x, dx, ds_tmp, p0_tmp, dv, xia0, Th, Ba, Ktt_tmp,
                        Ktt_tmp, ds_tmp, p0_tmp, M, C, K, R, G, p, Phi, J, &Vg,
                        &Kin);
  emxDestroyArray_real_T(J);
  emxDestroyArray_real_T(G);
  emxDestroyArray_real_T(R);
  emxDestroyArray_real_T(K);
  emxDestroyArray_real_T(C);
  emxDestroyArray_real_T(M);
  emxDestroyArray_real_T(Ba);
  emxDestroyArray_real_T(Th);
  emxDestroyArray_real_T(xia0);
  emxDestroyArray_real_T(dx);
  emxDestroyArray_real_T(x);
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_computeLagrangianFast();
  /* Terminate the application.
You do not need to do this more than one time. */
  computeLagrangianFast_terminate();
  return 0;
}

/* End of code generation (main.c) */
