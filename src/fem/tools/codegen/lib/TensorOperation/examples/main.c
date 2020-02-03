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
#include "rt_nonfinite.h"
#include "TensorOperation.h"
#include "main.h"
#include "TensorOperation_terminate.h"
#include "TensorOperation_initialize.h"

/* Function Declarations */
static void argInit_3x3_real_T(double result[9]);
static boolean_T argInit_boolean_T(void);
static double argInit_real_T(void);
static void main_TensorOperation(void);

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

static boolean_T argInit_boolean_T(void)
{
  return false;
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_TensorOperation(void)
{
  double dv0[9];
  double dv1[9];
  double T[36];

  /* Initialize function 'TensorOperation' input arguments. */
  /* Initialize function input argument 'A'. */
  /* Initialize function input argument 'B'. */
  /* Call the entry-point 'TensorOperation'. */
  argInit_3x3_real_T(dv0);
  argInit_3x3_real_T(dv1);
  TensorOperation(dv0, dv1, argInit_boolean_T(), T);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  TensorOperation_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_TensorOperation();

  /* Terminate the application.
     You do not need to do this more than one time. */
  TensorOperation_terminate();
  return 0;
}

/* End of code generation (main.c) */
