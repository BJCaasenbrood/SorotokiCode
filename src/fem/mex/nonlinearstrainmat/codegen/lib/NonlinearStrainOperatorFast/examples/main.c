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
#include "NonlinearStrainOperatorFast.h"
#include "NonlinearStrainOperatorFast_emxAPI.h"
#include "NonlinearStrainOperatorFast_terminate.h"
#include "NonlinearStrainOperatorFast_types.h"

/* Function Declarations */
static void argInit_3x3_real_T(double result[9]);

static emxArray_real_T *argInit_Unboundedx1_real_T(void);

static double argInit_real_T(void);

static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void);

static void c_main_NonlinearStrainOperatorF(void);

/* Function Definitions */
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

static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void)
{
  emxArray_real_T *result;
  int idx0;
  int idx1;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result = emxCreate_real_T(2, 2);
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }
  return result;
}

static void c_main_NonlinearStrainOperatorF(void)
{
  emxArray_real_T *Bg;
  emxArray_real_T *Bn;
  emxArray_real_T *N;
  emxArray_real_T *NN;
  emxArray_real_T *dNdx;
  double dv[9];
  double tau;
  emxInitArray_real_T(&Bn, 2);
  emxInitArray_real_T(&Bg, 2);
  emxInitArray_real_T(&NN, 2);
  /* Initialize function 'NonlinearStrainOperatorFast' input arguments. */
  /* Initialize function input argument 'N'. */
  N = argInit_Unboundedx1_real_T();
  /* Initialize function input argument 'dNdx'. */
  dNdx = c_argInit_UnboundedxUnbounded_r();
  /* Initialize function input argument 'F'. */
  /* Call the entry-point 'NonlinearStrainOperatorFast'. */
  argInit_3x3_real_T(dv);
  NonlinearStrainOperatorFast(N, dNdx, dv, Bn, Bg, NN, &tau);
  emxDestroyArray_real_T(NN);
  emxDestroyArray_real_T(Bg);
  emxDestroyArray_real_T(Bn);
  emxDestroyArray_real_T(dNdx);
  emxDestroyArray_real_T(N);
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  c_main_NonlinearStrainOperatorF();
  /* Terminate the application.
You do not need to do this more than one time. */
  NonlinearStrainOperatorFast_terminate();
  return 0;
}

/* End of code generation (main.c) */
