// VChooseK.c
// VCHOOSEK - Combinations of K elements [MEX]
// VCHOOSEK(V, K) creates a matrix, which rows are all combinations of
// choosing K elements of the vector V without order and without repititions.
//
// INPUT:
//   V: Array of class DOUBLE, SINGLE, (U)INT8/16/32/64, LOGICAL, CHAR.
//      The shape of the array does not matter.
//   K: Number of elements to choose. Scalar double with integer number.
//
// OUTPUT:
//   Y: Matrix of size [N!/K!(N-K)!, K] with N is the number of elements of V.
//      Y has the same class as the input V. If V has less than K elements, a
//      warning appears and the empty matrix is replied. As usual in Matlab,
//      empty inputs yield to empty Y.
//      A warning appears, if the output would exceed 500MB, and an error for
//      1GB. Both limits can be adjusted according to the available RAM in the
//      C-Mex source.
//
// NOTE: The output equals Matlab's NCHOOSEK, except that NCHOOSEK replies the
// number of combinations for a scalar V:
//   VChooseK(-1, 1) replies [-1]: One element taken out of a set of length one.
//   NCHOOSEK(-1, 1) fails at calculating N!/K!(N-K)!.
//
// EXAMPLES:
//   Choose 2 elements from [1,2,3,4]:
//     VChooseK(1:4, 2)  % ==> [1,2; 1,3; 1,4; 2,3; 2,4; 3,4]
//   For speed cast the input to integer types if possible:
//     Y = double(VChooseK(int16(1:1000), 2);
//   is faster than:
//     Y = VChooseK(1:1000, 2);
//   To get the combinations of cell arrays, use the combinations of the index:
//     C  = {'a', 'b', 'c', 'd'};
//     C2 = C(VChooseK(1:4, 2))
//     ==> C2 = {'a', 'b'; 'a', 'c'; 'a', 'd'; 'b', 'c'; 'b', 'd'; 'c', 'd'}
//
// COMPILE:
//   mex -O VChooseK.c
// Compatibility to 64-bit machines is assumed, but cannot be tested currently.
// On Linux the C99 comments must be considered (thanks Sebastiaan Breedveld):
//   mex CFLAGS="\$CFLAGS -std=C99" -O VChooseK.c
// Please run the unit-test TestVChooseK after compiling!
//
// INSPIRATION:
// Jos van der Geest has published an efficient NCHOOSE2, which is much faster
// than Matlab's NCHOOSEK:
//   http://www.mathworks.com/matlabcentral/fileexchange/20144
// I was curious, if a MEX version is much faster. At first I've created
// NChoose2q, which is 5 times faster than Jos' Matlab implementation. But the
// algorithm was so simple in C and did not consume temporary memory, that I
// expanded it to K=3,4,5 soon. The algorithm for free K was more challenging,
// but it needs just 3*K pointers for temporary memory. Therefore it is much
// faster than Matlab's NCHOOSEK (Matlab 6.5, 1.5GHz Pentium-M, 512 MB):
//                    NCHOOSEK:   VCHOOSEK:
//   (1:50, 5):       45.30 sec   0.32 sec      => 141 times faster
//   (1:20, 16):       0.52 sec   0.0024 sec    => 216 times faster
//   (int8(1:40), 4):  0.64 sec   0.001562 sec  => 410 times faster
//
// Related FEX related publications:
// (http://www.mathworks.com/matlabcentral/fileexchange/<number>)
//   NCTWO (Simone Scaringi): 20110
//   COMBINATOR (Matt Fig) [I prefer this for general tasks!]: 24325
//   PICK (Stefan Stoll) calls NCHOOSEK for unordered combinations: 12724
//   COMBINATIONS (Gautam Vallabha): 23080
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
//         Compilers: BCC5.5, LCC2.4/3.8, Open Watcom 1.8
//         Watcom lib is 20% faster than BCC and LCC libs for INT8!
// Author: Jan Simon, Heidelberg, (C) 2009 matlab.THIS_YEAR(a)nMINUSsimonDOTde
// License: BSD (use/copy/modify on own risk, but mention author)

/*
% $JRev: R0j V:007 Sum:5IALpldTXov9 Date:24-Dec-2008 02:59:35 $
% $File: User\JSim\Published\VChooseK\VChooseK.c $
% History:
% 001: 17-Dec-2009 08:10, MEX function for NCHOOSEK(X, 2).
% 007: 23-Dec-2008 17:40, Stop if output size exceed some limits.
%      Calculate output size as mwSize and DOUBLE to catch overflows.
*/

#include "mex.h"
#include <math.h>

// 32 bit array dimensions for Matlab 6.5:
// See option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef mwSize
#define mwSize  int
#define mwIndex int
#define MWSIZE_MAX 2147483647UL
#endif

// Limits for warning and error:
#define WARN_LIMIT   500000000L   // 500 MB, *must* be smaller than ERROR_LIMIT
#define ERROR_LIMIT 1000000000L   // 1 GB (recommended for 32&64 bit systems)

// Prototypes:
mwSize GetOutputLen(mwSize n, mwSize k, double *C_double);
void BadInputTypeError(void);

// The functions for each number of elements are absolutely equal except for
// the type of the array. But using some macros would result in ugly source
// code also. Using pointers to vectors of variable size would work, but is
// remarkably slower. Sorry for the pile of redundant code!
void Elem2_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem2_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem2_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem2_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void Elem3_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem3_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem3_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem3_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void Elem4_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem4_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem4_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem4_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void Elem5_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem5_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem5_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem5_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void ElemK_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY, mwSize K);
void ElemK_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY, mwSize K);
void ElemK_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY, mwSize K);
void ElemK_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY, mwSize K);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize nX, nY, K;
  const mxArray *X;
  double nY_double;
  int ElementSize;
  
  // Check number of inputs and outputs:
  if (nrhs != 2) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput", "2 inputs required.");
  }
  if (nlhs > 1) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput", "1 output allowed.");
  }
  
  // Check type of input:
  X = prhs[0];
  if (!mxIsNumeric(X) && !mxIsChar(X) & !mxIsLogical(X)) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput",
                       "Input array must be numerical, CHAR or LOGICAL.");
  }
  if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput",
                       "Input K must be a numerical scalar.");
  }
  
  // Get number of elements to choose and catch pathological cases:
  nX = mxGetNumberOfElements(X);
  K  = (mwSize) floor(mxGetScalar(prhs[1]) + 0.5);  // Rounding
  
  // Empty input => empty output:
  if (nX == 0) {
     plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     return;
  }
  
  // K == 1 => reply input as column vector, K == 0 => empty output:
  if (K <= 1) {
     if (K == 1) {  // Choose 1 element => output equals input:
        plhs[0] = mxDuplicateArray(X);
        mxSetM(plhs[0], nX);
        mxSetN(plhs[0], 1);
     } else if (K <= 0) {  // Empty matrix:
       plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     }
     return;
  }
  
  // Catch trivial K == nX and pathological K > nX:
  if (K >= nX) {
     if (K == nX) {
        plhs[0] = mxDuplicateArray(X);
        mxSetM(plhs[0], 1);
        mxSetN(plhs[0], nX);
     } else {
        mexWarnMsgIdAndTxt("JSimon:VChooseK:ShortInput",
                           "Too short input: K > N.");
        plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     }
     return;
  }
  
  // Calculate size of output (N!/(K!*(N-K)!)) and check limits - reduce
  // the risk of overflows, nY is calculated as DOUBLE additionally:
  nY          = GetOutputLen(nX, K, &nY_double);
  ElementSize = mxGetElementSize(X);
  nY_double  *= ElementSize;
  if (nY_double > WARN_LIMIT) {
     if (nY_double < ERROR_LIMIT && nY_double < MWSIZE_MAX) {
        mexWarnMsgIdAndTxt("JSimon:VChooseK:LargeOutput",
                           "Output will be large and take a long time.");
     } else {
        mexErrMsgIdAndTxt("JSimon:VChooseK:TooLargeOutput",
                          "Output would be too large.");
     }
  }
  
  // Create the output:
  // In a MEX file this stops with an error message automatically on failure.
  // For usage in an engine file, check result for NULL!
  plhs[0] = mxCreateNumericMatrix(nY, K, mxGetClassID(X), mxREAL);
  
  // Call different core functions according to the element size:
  switch (K) {
    case 2:  // K == 2:
      switch (ElementSize) {
        case 8:  Elem2_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem2_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem2_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem2_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;
        
    case 3:  // K == 3:
      switch (ElementSize) {
        case 8:  Elem3_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem3_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem3_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem3_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;
        
    case 4:  // K == 4:
      switch (ElementSize) {
        case 8:  Elem4_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem4_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem4_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem4_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;
      
    case 5:  // K == 5:
      switch (ElementSize) {
        case 8:  Elem5_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem5_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem5_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem5_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;

    default:  // K > 5
      switch (ElementSize) {
        case 8:  ElemK_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY, K);
                 break;
        case 4:  ElemK_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY, K);
                 break;
        case 2:  ElemK_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY, K);
                 break;
        case 1:  ElemK_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY, K);
                 break;
        default: BadInputTypeError();
      }
  }
  
  return;
}

// Error about bad input type: =================================================
void BadInputTypeError(void)
{
  mexErrMsgIdAndTxt("JSimon:VChooseK:BadInputType",
                    "Input must have 1, 2, 4 or 8 bytes per element.");
}

// Get output size: ============================================================
mwSize GetOutputLen(mwSize n, mwSize k, double *C_double)
{
  // Number of samples for taking k elements from a set of n.
  // Calculate the exact value as INT32 (or INT64 on 64 bit machine) and the
  // eventually less precise DOUBLE value, which is less susceptible for an
  // overflow.
  mwSize i, ai, bi;
  double ad, bd;
  
  bi = n - k;
  ai = bi + 1;
  ad = (double) ai;
  bd = (double) bi;
  for (i = 2; i <= k; i++) {
     ai += (ai * bi) / i;
     ad += (ad * bd) / i;
  }
  
  *C_double = ad;
  
  return (ai);
}

// 2 elements: =================================================================
void Elem2_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  // Called for input data with 8 bytes per element: DOUBLE, (U)INT64.
  // 2 elements - fill 2 columns of output simultaneously.
  double *Xq, *Xf, *Y1, *Y2;
  
  Y1 = Yp;             // 1st column of output
  Y2 = Y1 + nY;        // 2nd column of output
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {  // With repititions: Xq = Xp;
      *Y1++ = *Xp;     // 1st column: slow changing
      *Y2++ = *Xq;     // 2nd column: fast changing
    }
  }

  return;
}

// =============================================================================
void Elem2_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  // Called for input data with 4 bytes per element: SINGLE, (U)INT32.
  // 2 elements - fill 2 columns of output simultaneously.
  INT32_T *Xq, *Xf, *Y1, *Y2;

  Y1 = Yp;             // 1st column of output
  Y2 = Y1 + nY;        // 2nd column of output
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {
      *Y1++ = *Xp;     // 1st column: slow changing
      *Y2++ = *Xq;     // 2nd column: fast changing
    }
  }

  return;
}

// =============================================================================
void Elem2_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  // Called for input data with 2 bytes per element: (U)INT16, CHAR.
  // 2 elements - fill 2 columns of output simultaneously.
  INT16_T *Xq, *Xf, *Y1, *Y2;

  Y1 = Yp;             // 1st column of output
  Y2 = Y1 + nY;        // 2nd column of output
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {
      *Y1++ = *Xp;     // 1st column: slow changing
      *Y2++ = *Xq;     // 2nd column: fast changing
    }
  }

  return;
}

// =============================================================================
void Elem2_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  // Called for input data with 1 byte per element: (U)INT8, LOGICAL.
  // 2 elements - fill 2 columns of output simultaneously.
  INT8_T *Xq, *Xf, *Y1, *Y2;

  Y1 = Yp;             // 1st column of output
  Y2 = Y1 + nY;        // 2nd column of output
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {
      *Y1++ = *Xp;     // 1st column: slow changing
      *Y2++ = *Xq;     // 2nd column: fast changing
    }
  }

  return;
}

// 3 elements: =================================================================
void Elem3_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  // Called for input data with 8 bytes per element: DOUBLE, (U)INT64.
  // Fill 3 columns of output simultaneously.
  double *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           // 1st column of output
  Y2  = Y1 + nY;      // 2nd column of output
  Y3  = Y2 + nY;      // 3rd column of output
  X3f = Xp + nX;      // 3rd column ends with element nX
  X2f = X3f - 1;      // 2nd column ends with element nX - 1
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {          // "<", not "<=" !
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  // 1st column: slow changing
        *Y2++ = *Xq;
        *Y3++ = *Xr;  // 3rd column: fastest changing
      }
    }
  }

  return;
}

// =============================================================================
void Elem3_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  // Called for input data with 4 bytes per element: SINGLE, (U)INT32.
  // Fill 3 columns of output simultaneously.
  INT32_T *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           // 1st column of output
  Y2  = Y1 + nY;      // 2nd column of output
  Y3  = Y2 + nY;      // 3rd column of output
  X3f = Xp + nX;      // 3rd column ends with element nX
  X2f = X3f - 1;      // 2nd column ends with element nX - 1
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {          // "<", not "<=" !
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  // 1st column: slow changing
        *Y2++ = *Xq;
        *Y3++ = *Xr;  // 3rd column: fastest changing
      }
    }
  }
  
  return;
}

// =============================================================================
void Elem3_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  // Called for input data with 2 bytes per element: (U)INT16, CHAR.
  // Fill 3 columns of output simultaneously.
  INT16_T *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           // 1st column of output
  Y2  = Y1 + nY;      // 2nd column of output
  Y3  = Y2 + nY;      // 3rd column of output
  X3f = Xp + nX;      // 3rd column ends with element nX
  X2f = X3f - 1;      // 2nd column ends with element nX - 1
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {          // "<", not "<=" !
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  // 1st column: slow changing
        *Y2++ = *Xq;
        *Y3++ = *Xr;  // 3rd column: fastest changing
      }
    }
  }
  
  return;
}

// =============================================================================
void Elem3_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  // Called for input data with 1 byte per element: (U)INT8, LOGICAL.
  // Fill 3 columns of output simultaneously.
  INT8_T *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           // 1st column of output
  Y2  = Y1 + nY;      // 2nd column of output
  Y3  = Y2 + nY;      // 3rd column of output
  X3f = Xp + nX;      // 3rd column ends with element nX
  X2f = X3f - 1;      // 2nd column ends with element nX - 1
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {          // "<", not "<=" !
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  // 1st column: slow changing
        *Y2++ = *Xq;
        *Y3++ = *Xr;  // 3rd column: fastest changing
      }
    }
  }
  
  return;
}

// 4 elements: =================================================================
void Elem4_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  // Called for input data with 8 bytes per element: DOUBLE, (U)INT64.
  // Fill 4 columns of output simultaneously.
  double *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;              // 1st column of output
  Y2  = Y1 + nY;         // 2nd column
  Y3  = Y2 + nY;         // 3rd column
  Y4  = Y3 + nY;         // 4th column
  X4f = Xp + nX;         // 4th column ends with element nX
  X2f = X4f - 2;         // 2nd column ends with element nX - 2
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {       // "<", not "<="
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;  // 1st column: slowest changing
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs;  // 4th column: fastest changing
        }
      }
    }
  }

  return;
}

// =============================================================================
void Elem4_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  // Called for input data with 4 bytes per element: SINGLE, (U)INT32.
  // Fill 4 columns of output simultaneously.
  INT32_T *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;              // 1st column of output
  Y2  = Y1 + nY;         // 2nd column
  Y3  = Y2 + nY;         // 3rd column
  Y4  = Y3 + nY;         // 4th column
  X4f = Xp + nX;         // 4th column ends with element nX
  X2f = X4f - 2;         // 2nd column ends with element nX - 2
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {       // "<", not "<="
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;  // 1st column: slowest changing
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs;  // 4th column: fastest changing
        }
      }
    }
  }
  
  return;
}

// =============================================================================
void Elem4_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  // Called for input data with 2 bytes per element: (U)INT16, CHAR.
  // Fill 4 columns of output simultaneously.
  INT16_T *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;              // 1st column of output
  Y2  = Y1 + nY;         // 2nd column
  Y3  = Y2 + nY;         // 3rd column
  Y4  = Y3 + nY;         // 4th column
  X4f = Xp + nX;         // 4th column ends with element nX
  X2f = X4f - 2;         // 2nd column ends with element nX - 2
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {       // "<", not "<="
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;  // 1st column: slowest changing
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs;  // 4th column: fastest changing
        }
      }
    }
  }
  
  return;
}

// =============================================================================
void Elem4_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  // Called for input data with 1 byte per element: (U)INT8, LOGICAL.
  // Fill 4 columns of output simultaneously.
  INT8_T *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;              // 1st column of output
  Y2  = Y1 + nY;         // 2nd column
  Y3  = Y2 + nY;         // 3rd column
  Y4  = Y3 + nY;         // 4th column
  X4f = Xp + nX;         // 4th column ends with element nX
  X2f = X4f - 2;         // 2nd column ends with element nX - 2
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {       // "<", not "<="
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;  // 1st column: slowest changing
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs;  // 4th column: fastest changing
        }
      }
    }
  }
  
  return;
}

// 5 elements: =================================================================
void Elem5_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  // Called for input data with 8 bytes per element: DOUBLE, (U)INT64.
  // Fill 5 columns of output simultaneously.
  double *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
         *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               // 1st column of output
  Y2  = Y1 + nY;          // 2nd column
  Y3  = Y2 + nY;          // 3rd column
  Y4  = Y3 + nY;          // 4th column
  Y5  = Y4 + nY;          // 5th column
  X5f = Xp + nX;          // 5th column ends with element nX
                          // Used for 4th column also
  X3f = X5f - 2;          // 3rd column ends with element nX - 2
                          // Used for 2nd column also
  X1f = X3f - 2;          // 1st column ends with element nX - 4
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp;  // 1st column: slowest changing
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  // 5th column: fastest changing
          }
        }
      }
    }
  }

  return;
}

// =============================================================================
void Elem5_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  // Called for input data with 4 bytes per element: SINGLE, (U)INT32.
  // Fill 5 columns of output simultaneously.
  INT32_T *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
          *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               // 1st column of output
  Y2  = Y1 + nY;          // 2nd column
  Y3  = Y2 + nY;          // 3rd column
  Y4  = Y3 + nY;          // 4th column
  Y5  = Y4 + nY;          // 5th column
  X5f = Xp + nX;          // 5th column ends with element nX
  X3f = X5f - 2;          // 3rd column ends with element nX - 2
  X1f = X3f - 2;          // 1st column ends with element nX - 4
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp;  // 1st column: slowest changing
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  // 5th column: fastest changing
          }
        }
      }
    }
  }

  return;
}

// =============================================================================
void Elem5_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  // Called for input data with 2 bytes per element: (U)INT16, CHAR.
  // Fill 5 columns of output simultaneously.
  INT16_T *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
          *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               // 1st column of output
  Y2  = Y1 + nY;          // 2nd column
  Y3  = Y2 + nY;          // 3rd column
  Y4  = Y3 + nY;          // 4th column
  Y5  = Y4 + nY;          // 5th column
  X5f = Xp + nX;          // 5th column ends with element nX
  X3f = X5f - 2;          // 3rd column ends with element nX - 2
  X1f = X3f - 2;          // 1st column ends with element nX - 4
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp;  // 1st column: slowest changing
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  // 5th column: fastest changing
          }
        }
      }
    }
  }
  
  return;
}

// =============================================================================
void Elem5_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  // Called for input data with 1 byte per element: (U)INT8, LOGICAL.
  // Fill 5 columns of output simultaneously.
  INT8_T *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
         *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               // 1st column of output
  Y2  = Y1 + nY;          // 2nd column
  Y3  = Y2 + nY;          // 3rd column
  Y4  = Y3 + nY;          // 4th column
  Y5  = Y4 + nY;          // 5th column
  X5f = Xp + nX;          // 5th column ends with element nX
  X3f = X5f - 2;          // 3rd column ends with element nX - 2
  X1f = X3f - 2;          // 1st column ends with element nX - 4
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp;  // 1st column: slowest changing
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  // 5th column: fastest changing
          }
        }
      }
    }
  }
  
  return;
}

// =============================================================================
// Of course we could copy and paste the Elem[2/3/4/5]_NByte code until we reach
// the outer limits, but this is boring. So let's create a general algorithm,
// although it is clear, that such algorithms must have been implemented
// anywhere already?!
// =============================================================================
void ElemK_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY, mwSize K)
{
  // Called for input data with 8 bytes per element: DOUBLE, UINT64.
  // Fill all columns of output simultaneously.
  double *Xq, *Xr,
         **A, **Ap, **Af, *Yq,
         **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
    
  // Array of input pointers [B] and final indices for each column [F]:
  if ((B = (double **) mxMalloc(K * sizeof(double *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (double **) mxMalloc(K * sizeof(double *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  // Set K input pointers and final index for each column:
  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  // Array of output pointers:
  if ((A = (double **) mxMalloc(K * sizeof(double *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  // Set pointers to the K output columns of Y:
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  // Now the loops of the standard method (e.g. Elem3_8Byte) are simulated:
  // The pointers formerly contained in a cascade of loops, are stored in the
  // array B. The final index for each "loop" is stored in F.
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    // Copy values from the input to left K-1 columns, just the value of the
    // K.th column changes:
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       // Copy values from input to output for all columns:
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;           // Combine this to "*(*Ap)++ = **Bp++;" ?!
       }
       (*Bq)++;             // Advance input pointer for K.th column only
    }
    
    // Find the right-most not exhausted loop:
    do {
       Bq--;
       Fp--;
       if (Bq < B) {        // Even the 1st column is exhausted - ready!
          mxFree(A);        // Thanks for the fish
          mxFree(B);
          mxFree(F);
          return;           // Goodbye
       }
    } while (*Bq >= *Fp);   // Input pointer is before end of loop
    
    (*Bq)++;                // Advance pointer of last not exhausted loop
    while (++Bq < Bf) {     // Initialize the following loops:
      *Bq = *(Bq - 1) + 1;  // With repititions: *Bq = *(Bq - 1)
    }
  }  // while(1) loop
  
  // Never reached this line!
}

// =============================================================================
void ElemK_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY, mwSize K)
{
  // Called for input data with 4 bytes per element: (U)INT32, SINGLE.
  // Fill all columns of output simultaneously.
  INT32_T *Xq, *Xr,
          **A, **Ap, **Af, *Yq,
          **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
    
  // Array of input pointers [B] and final indices for each column [F]:
  if ((B = (INT32_T **) mxMalloc(K * sizeof(INT32_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (INT32_T **) mxMalloc(K * sizeof(INT32_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  // Set K input pointers and final index for each column:
  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  // Array of output pointers:
  if ((A = (INT32_T **) mxMalloc(K * sizeof(INT32_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  // Set pointers to the K output columns of Y:
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  // Now the loops of the standard method (e.g. Elem3_8Byte) are simulated:
  // The pointers formerly contained in a cascade of loops, are stored in the
  // array B. The final index for each "loop" is stored in F.
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    // Copy values from the input to left K-1 columns, just the value of the
    // K.th column changes:
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       // Copy values from input to output for all columns:
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;           // Combine this to "*(*Ap)++ = **Bp++;" ?!
       }
       (*Bq)++;             // Advance input pointer for K.th column only
    }
    
    // Find the right-most not exhausted loop:
    do {
       Bq--;
       Fp--;
       if (Bq < B) {        // Even the 1st column is exhausted - ready!
          mxFree(A);        // Thanks for the fish
          mxFree(B);
          mxFree(F);
          return;           // Goodbye
       }
    } while (*Bq >= *Fp);   // Input pointer is before end of loop
    
    (*Bq)++;                // Advance pointer of last not exhausted loop
    while (++Bq < Bf) {     // Initialize the following loops:
      *Bq = *(Bq - 1) + 1;  // With repititions: *Bq = *(Bq - 1)
    }
  }  // while(1) loop
  
  // Never reached this line!
}

// =============================================================================
void ElemK_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY, mwSize K)
{
  // Called for input data with 2 bytes per element: (U)INT16, CHAR.
  // Fill all columns of output simultaneously.
  INT16_T *Xq, *Xr,
          **A, **Ap, **Af, *Yq,
          **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
  
  // Array of input pointers [B] and final indices for each column [F]:
  if ((B = (INT16_T **) mxMalloc(K * sizeof(INT16_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (INT16_T **) mxMalloc(K * sizeof(INT16_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  // Set K input pointers and final index for each column:
  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  // Array of output pointers:
  if ((A = (INT16_T **) mxMalloc(K * sizeof(INT16_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  // Set pointers to the K output columns of Y:
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  // Now the loops of the standard method (e.g. Elem3_8Byte) are simulated:
  // The pointers formerly contained in a cascade of loops, are stored in the
  // array B. The final index for each "loop" is stored in F.
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    // Copy values from the input to left K-1 columns, just the value of the
    // K.th column changes:
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       // Copy values from input to output for all columns:
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;           // Combine this to "*(*Ap)++ = **Bp++;" ?!
       }
       (*Bq)++;             // Advance input pointer for K.th column only
    }
    
    // Find the right-most not exhausted loop:
    do {
       Bq--;
       Fp--;
       if (Bq < B) {        // Even the 1st column is exhausted - ready!
          mxFree(A);        // Thanks for the fish
          mxFree(B);
          mxFree(F);
          return;           // Goodbye
       }
    } while (*Bq >= *Fp);   // Input pointer is before end of loop
    
    (*Bq)++;                // Advance pointer of last not exhausted loop
    while (++Bq < Bf) {     // Initialize the following loops:
      *Bq = *(Bq - 1) + 1;  // With repititions: *Bq = *(Bq - 1)
    }
  }  // while(1) loop
  
  // Never reached this line!
}

// =============================================================================
void ElemK_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY, mwSize K)
{
  // Called for input data with 1 byte per element: (U)INT8, LOGICAL.
  // Fill all columns of output simultaneously.
  INT8_T  *Xq, *Xr,
          **A, **Ap, **Af, *Yq,
          **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
    
  // Array of input pointers [B] and final indices for each column [F]:
  if ((B = (INT8_T **) mxMalloc(K * sizeof(INT8_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (INT8_T **) mxMalloc(K * sizeof(INT8_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  // Set K input pointers and final index for each column:
  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  // Array of output pointers:
  if ((A = (INT8_T **) mxMalloc(K * sizeof(INT8_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  // Set pointers to the K output columns of Y:
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  // Now the loops of the standard method (e.g. Elem3_8Byte) are simulated:
  // The pointers formerly contained in a cascade of loops, are stored in the
  // array B. The final index for each "loop" is stored in F.
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    // Copy values from the input to left K-1 columns, just the value of the
    // K.th column changes:
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       // Copy values from input to output for all columns:
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;           // Combine this to "*(*Ap)++ = **Bp++;" ?!
       }
       (*Bq)++;             // Advance input pointer for K.th column only
    }
    
    // Find the right-most not exhausted loop:
    do {
       Bq--;
       Fp--;
       if (Bq < B) {        // Even the 1st column is exhausted - ready!
          mxFree(A);        // Thanks for the fish
          mxFree(B);
          mxFree(F);
          return;           // Goodbye
       }
    } while (*Bq >= *Fp);   // Input pointer is before end of loop
    
    (*Bq)++;                // Advance pointer of last not exhausted loop
    while (++Bq < Bf) {     // Initialize the following loops:
      *Bq = *(Bq - 1) + 1;  // With repititions: *Bq = *(Bq - 1)
    }
  }  // while(1) loop
  
  // Never reached this line!
}
