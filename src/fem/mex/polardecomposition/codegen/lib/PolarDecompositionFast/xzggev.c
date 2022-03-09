/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzggev.c
 *
 * Code generation for function 'xzggev'
 *
 */

/* Include files */
#include "xzggev.h"
#include "PolarDecompositionFast_rtwutil.h"
#include "rt_nonfinite.h"
#include "xzhgeqz.h"
#include "xzlartg.h"
#include "xztgevc.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void xzggev(creal_T A[9], int *info, creal_T alpha1[3], creal_T beta1[3],
            creal_T V[9])
{
  creal_T atmp;
  double absxk;
  double anrm;
  double anrmto;
  double c;
  double cfrom1;
  double cto1;
  double ctoc;
  double d;
  double d1;
  double d2;
  double d3;
  int rscale[3];
  int A_tmp;
  int exitg2;
  int exitg3;
  int i;
  int ihi;
  int ii;
  int ilo;
  int j;
  int k;
  int nzcount;
  signed char b_I[9];
  boolean_T exitg1;
  boolean_T exitg4;
  boolean_T guard1 = false;
  boolean_T ilascl;
  boolean_T notdone;
  *info = 0;
  anrm = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 9)) {
    absxk = rt_hypotd_snf(A[k].re, A[k].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }
      k++;
    }
  }
  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    alpha1[0].re = rtNaN;
    alpha1[0].im = 0.0;
    beta1[0].re = rtNaN;
    beta1[0].im = 0.0;
    alpha1[1].re = rtNaN;
    alpha1[1].im = 0.0;
    beta1[1].re = rtNaN;
    beta1[1].im = 0.0;
    alpha1[2].re = rtNaN;
    alpha1[2].im = 0.0;
    beta1[2].re = rtNaN;
    beta1[2].im = 0.0;
    for (A_tmp = 0; A_tmp < 9; A_tmp++) {
      V[A_tmp].re = rtNaN;
      V[A_tmp].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
      guard1 = true;
    } else if (anrm > 1.4885657073574029E+138) {
      anrmto = 1.4885657073574029E+138;
      ilascl = true;
      guard1 = true;
    }
    if (guard1) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          c = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          c = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          c = ctoc / absxk;
          notdone = false;
        }
        for (A_tmp = 0; A_tmp < 9; A_tmp++) {
          A[A_tmp].re *= c;
          A[A_tmp].im *= c;
        }
      }
    }
    rscale[0] = 1;
    rscale[1] = 1;
    rscale[2] = 1;
    ilo = 1;
    ihi = 3;
    do {
      exitg3 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg1 = false;
      while ((!exitg1) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        k = 0;
        exitg4 = false;
        while ((!exitg4) && (k <= ihi - 1)) {
          A_tmp = (ii + 3 * k) - 1;
          if ((A[A_tmp].re != 0.0) || (A[A_tmp].im != 0.0) || (ii == k + 1)) {
            if (nzcount == 0) {
              j = k + 1;
              nzcount = 1;
              k++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            k++;
          }
        }
        if (nzcount < 2) {
          notdone = true;
          exitg1 = true;
        } else {
          ii--;
        }
      }
      if (!notdone) {
        exitg3 = 2;
      } else {
        if (i != ihi) {
          atmp = A[i - 1];
          A[i - 1] = A[ihi - 1];
          A[ihi - 1] = atmp;
          atmp = A[i + 2];
          A[i + 2] = A[ihi + 2];
          A[ihi + 2] = atmp;
          atmp = A[i + 5];
          A[i + 5] = A[ihi + 5];
          A[ihi + 5] = atmp;
        }
        if (j != ihi) {
          for (k = 0; k < ihi; k++) {
            ii = k + 3 * (j - 1);
            atmp = A[ii];
            A_tmp = k + 3 * (ihi - 1);
            A[ii] = A[A_tmp];
            A[A_tmp] = atmp;
          }
        }
        rscale[ihi - 1] = j;
        ihi--;
        if (ihi == 1) {
          rscale[0] = 1;
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);
    if (exitg3 != 1) {
      do {
        exitg2 = 0;
        i = 0;
        j = 0;
        notdone = false;
        k = ilo;
        exitg1 = false;
        while ((!exitg1) && (k <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = k;
          ii = ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= ihi)) {
            A_tmp = (ii + 3 * (k - 1)) - 1;
            if ((A[A_tmp].re != 0.0) || (A[A_tmp].im != 0.0) || (ii == k)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                ii++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              ii++;
            }
          }
          if (nzcount < 2) {
            notdone = true;
            exitg1 = true;
          } else {
            k++;
          }
        }
        if (!notdone) {
          exitg2 = 1;
        } else {
          if (i != ilo) {
            for (k = ilo; k < 4; k++) {
              ii = 3 * (k - 1);
              nzcount = (i + ii) - 1;
              atmp = A[nzcount];
              A_tmp = (ilo + ii) - 1;
              A[nzcount] = A[A_tmp];
              A[A_tmp] = atmp;
            }
          }
          if (j != ilo) {
            for (k = 0; k < ihi; k++) {
              ii = k + 3 * (j - 1);
              atmp = A[ii];
              A_tmp = k + 3 * (ilo - 1);
              A[ii] = A[A_tmp];
              A[A_tmp] = atmp;
            }
          }
          rscale[ilo - 1] = j;
          ilo++;
          if (ilo == ihi) {
            rscale[ilo - 1] = ilo;
            exitg2 = 1;
          }
        }
      } while (exitg2 == 0);
    }
    for (A_tmp = 0; A_tmp < 9; A_tmp++) {
      b_I[A_tmp] = 0;
    }
    b_I[0] = 1;
    b_I[4] = 1;
    b_I[8] = 1;
    for (A_tmp = 0; A_tmp < 9; A_tmp++) {
      V[A_tmp].re = b_I[A_tmp];
      V[A_tmp].im = 0.0;
    }
    if (ihi >= ilo + 2) {
      ii = ilo;
      while (ii < 2) {
        xzlartg(A[1], A[2], &c, &atmp, &A[1]);
        A[2].re = 0.0;
        A[2].im = 0.0;
        cfrom1 = atmp.re * A[5].re - atmp.im * A[5].im;
        cto1 = atmp.re * A[5].im + atmp.im * A[5].re;
        absxk = c * A[4].re;
        ctoc = c * A[4].im;
        d = A[4].im;
        d1 = A[4].re;
        A[5].re = c * A[5].re - (atmp.re * A[4].re + atmp.im * A[4].im);
        A[5].im = c * A[5].im - (atmp.re * d - atmp.im * d1);
        A[4].re = absxk + cfrom1;
        A[4].im = ctoc + cto1;
        cfrom1 = atmp.re * A[8].re - atmp.im * A[8].im;
        cto1 = atmp.re * A[8].im + atmp.im * A[8].re;
        absxk = c * A[7].re;
        ctoc = c * A[7].im;
        d = A[7].im;
        d1 = A[7].re;
        A[8].re = c * A[8].re - (atmp.re * A[7].re + atmp.im * A[7].im);
        A[8].im = c * A[8].im - (atmp.re * d - atmp.im * d1);
        A[7].re = absxk + cfrom1;
        A[7].im = ctoc + cto1;
        atmp.re = -atmp.re;
        atmp.im = -atmp.im;
        absxk = c * A[6].re + (atmp.re * A[3].re - atmp.im * A[3].im);
        ctoc = c * A[6].im + (atmp.re * A[3].im + atmp.im * A[3].re);
        d = A[6].im;
        d1 = A[6].re;
        A[3].re = c * A[3].re - (atmp.re * A[6].re + atmp.im * A[6].im);
        A[3].im = c * A[3].im - (atmp.re * d - atmp.im * d1);
        A[6].re = absxk;
        A[6].im = ctoc;
        cfrom1 = atmp.re * A[4].re - atmp.im * A[4].im;
        cto1 = atmp.re * A[4].im + atmp.im * A[4].re;
        d2 = c * A[7].re;
        d3 = c * A[7].im;
        d = A[7].im;
        d1 = A[7].re;
        A[4].re = c * A[4].re - (atmp.re * A[7].re + atmp.im * A[7].im);
        A[4].im = c * A[4].im - (atmp.re * d - atmp.im * d1);
        A[7].re = d2 + cfrom1;
        A[7].im = d3 + cto1;
        cfrom1 = atmp.re * A[5].re - atmp.im * A[5].im;
        cto1 = atmp.re * A[5].im + atmp.im * A[5].re;
        d2 = c * A[8].re;
        d3 = c * A[8].im;
        d = A[8].im;
        d1 = A[8].re;
        A[5].re = c * A[5].re - (atmp.re * A[8].re + atmp.im * A[8].im);
        A[5].im = c * A[5].im - (atmp.re * d - atmp.im * d1);
        A[8].re = d2 + cfrom1;
        A[8].im = d3 + cto1;
        cfrom1 = atmp.re * V[3].re - atmp.im * V[3].im;
        cto1 = atmp.re * V[3].im + atmp.im * V[3].re;
        V[3].re = c * V[3].re - (atmp.re * V[6].re + atmp.im * V[6].im);
        V[3].im = c * V[3].im - (atmp.re * V[6].im - atmp.im * V[6].re);
        V[6].re = c * V[6].re + cfrom1;
        V[6].im = c * V[6].im + cto1;
        cfrom1 = atmp.re * V[4].re - atmp.im * V[4].im;
        cto1 = atmp.re * V[4].im + atmp.im * V[4].re;
        V[4].re = c * V[4].re - (atmp.re * V[7].re + atmp.im * V[7].im);
        V[4].im = c * V[4].im - (atmp.re * V[7].im - atmp.im * V[7].re);
        V[7].re = c * V[7].re + cfrom1;
        V[7].im = c * V[7].im + cto1;
        cfrom1 = atmp.re * V[5].re - atmp.im * V[5].im;
        cto1 = atmp.re * V[5].im + atmp.im * V[5].re;
        V[5].re = c * V[5].re - (atmp.re * V[8].re + atmp.im * V[8].im);
        V[5].im = c * V[5].im - (atmp.re * V[8].im - atmp.im * V[8].re);
        V[8].re = c * V[8].re + cfrom1;
        V[8].im = c * V[8].im + cto1;
        ii = 2;
      }
    }
    xzhgeqz(A, ilo, ihi, V, info, alpha1, beta1);
    if (*info == 0) {
      xztgevc(A, V);
      if (ilo > 1) {
        for (i = ilo - 2; i + 1 >= 1; i--) {
          k = rscale[i] - 1;
          if (rscale[i] != i + 1) {
            atmp = V[i];
            V[i] = V[k];
            V[k] = atmp;
            atmp = V[i + 3];
            V[i + 3] = V[k + 3];
            V[k + 3] = atmp;
            atmp = V[i + 6];
            V[i + 6] = V[k + 6];
            V[k + 6] = atmp;
          }
        }
      }
      if (ihi < 3) {
        A_tmp = ihi + 1;
        for (i = A_tmp; i < 4; i++) {
          ii = rscale[i - 1];
          if (ii != i) {
            atmp = V[i - 1];
            V[i - 1] = V[ii - 1];
            V[ii - 1] = atmp;
            atmp = V[i + 2];
            V[i + 2] = V[ii + 2];
            V[ii + 2] = atmp;
            atmp = V[i + 5];
            V[i + 5] = V[ii + 5];
            V[ii + 5] = atmp;
          }
        }
      }
      for (k = 0; k < 3; k++) {
        d = V[3 * k].re;
        d1 = V[3 * k].im;
        absxk = fabs(d) + fabs(d1);
        ii = 3 * k + 1;
        d2 = V[ii].re;
        d3 = V[ii].im;
        ctoc = fabs(d2) + fabs(d3);
        if (ctoc > absxk) {
          absxk = ctoc;
        }
        nzcount = 3 * k + 2;
        cfrom1 = V[nzcount].re;
        cto1 = V[nzcount].im;
        ctoc = fabs(cfrom1) + fabs(cto1);
        if (ctoc > absxk) {
          absxk = ctoc;
        }
        if (absxk >= 6.7178761075670888E-139) {
          absxk = 1.0 / absxk;
          d *= absxk;
          V[3 * k].re = d;
          d1 *= absxk;
          V[3 * k].im = d1;
          d2 *= absxk;
          V[ii].re = d2;
          d3 *= absxk;
          V[ii].im = d3;
          cfrom1 *= absxk;
          V[nzcount].re = cfrom1;
          cto1 *= absxk;
          V[nzcount].im = cto1;
        }
      }
      if (ilascl) {
        notdone = true;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            c = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            c = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            c = anrm / anrmto;
            notdone = false;
          }
          alpha1[0].re *= c;
          alpha1[0].im *= c;
          alpha1[1].re *= c;
          alpha1[1].im *= c;
          alpha1[2].re *= c;
          alpha1[2].im *= c;
        }
      }
    }
  }
}

/* End of code generation (xzggev.c) */
