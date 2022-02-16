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
#include "computeForwardKinematicsFast.h"
#include "computeForwardKinematicsFast_emxAPI.h"
#include "computeForwardKinematicsFast_terminate.h"
#include "computeForwardKinematicsFast_types.h"

/* Function Declarations */
static void argInit_3x1_real_T(double result[3]);

static void argInit_3x3_real_T(double result[9]);

static emxArray_real_T *argInit_6x1xUnbounded_real_T(void);

static emxArray_real_T *argInit_6xUnbounded_real_T(void);

static emxArray_real_T *argInit_Unboundedx1_real_T(void);

static double argInit_real_T(void);

static emxArray_real_T *c_argInit_UnboundedxUnboundedxU(void);

static void c_main_computeForwardKinematics(void);

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

static void c_main_computeForwardKinematics(void)
{
  emxArray_real_T *Ba;
  emxArray_real_T *Jtmp;
  emxArray_real_T *Th;
  emxArray_real_T *gtmp;
  emxArray_real_T *x;
  emxArray_real_T *xia0;
  double dv1[9];
  double dv[3];
  emxInitArray_real_T(&gtmp, 3);
  emxInitArray_real_T(&Jtmp, 3);
  /* Initialize function 'computeForwardKinematicsFast' input arguments. */
  /* Initialize function input argument 'x'. */
  x = argInit_Unboundedx1_real_T();
  /* Initialize function input argument 'p0'. */
  /* Initialize function input argument 'Phi0'. */
  /* Initialize function input argument 'xia0'. */
  xia0 = argInit_6x1xUnbounded_real_T();
  /* Initialize function input argument 'Th'. */
  Th = c_argInit_UnboundedxUnboundedxU();
  /* Initialize function input argument 'Ba'. */
  Ba = argInit_6xUnbounded_real_T();
  /* Call the entry-point 'computeForwardKinematicsFast'. */
  argInit_3x1_real_T(dv);
  argInit_3x3_real_T(dv1);
  computeForwardKinematicsFast(x, argInit_real_T(), dv, dv1, xia0, Th, Ba, gtmp,
                               Jtmp);
  emxDestroyArray_real_T(Jtmp);
  emxDestroyArray_real_T(gtmp);
  emxDestroyArray_real_T(Ba);
  emxDestroyArray_real_T(Th);
  emxDestroyArray_real_T(xia0);
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
  c_main_computeForwardKinematics();
  /* Terminate the application.
You do not need to do this more than one time. */
  computeForwardKinematicsFast_terminate();
  return 0;
}

/* End of code generation (main.c) */
