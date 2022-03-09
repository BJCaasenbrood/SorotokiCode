/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xhseqr.c
 *
 * Code generation for function 'xhseqr'
 *
 */

/* Include files */
#include "xhseqr.h"
#include "rt_nonfinite.h"
#include "xdlanv2.h"
#include "xzlarfg.h"
#include <math.h>

/* Function Definitions */
int xhseqr(double h[9], double z[9])
{
  double v[3];
  double aa;
  double ab;
  double ba;
  double bb;
  double d;
  double h22;
  double rt1r;
  double s;
  double tst;
  int b_k;
  int hoffset;
  int i;
  int info;
  int its;
  int j;
  int k;
  int m;
  int nr;
  int sum1_tmp_tmp;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T goto150;
  info = 0;
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
  h[2] = 0.0;
  i = 2;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    k = -1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = i - 1;
      exitg3 = false;
      while ((!exitg3) && (k + 2 > 1)) {
        sum1_tmp_tmp = k + 3 * k;
        ba = fabs(h[sum1_tmp_tmp + 1]);
        if (ba <= 3.0062525400134592E-292) {
          exitg3 = true;
        } else {
          nr = 3 * (k + 1);
          hoffset = k + nr;
          bb = fabs(h[hoffset + 1]);
          tst = fabs(h[sum1_tmp_tmp]) + bb;
          if (tst == 0.0) {
            if (k >= 1) {
              tst = fabs(h[k + 3 * (k - 1)]);
            }
            if (k + 3 <= 3) {
              tst += fabs(h[nr + 2]);
            }
          }
          if (ba <= 2.2204460492503131E-16 * tst) {
            tst = fabs(h[hoffset]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }
            tst = fabs(h[sum1_tmp_tmp] - h[hoffset + 1]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }
            s = aa + ab;
            if (ba * (ab / s) <=
                fmax(3.0062525400134592E-292,
                     2.2204460492503131E-16 * (bb * (aa / s)))) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }
      if (k + 2 > 1) {
        h[(k + 3 * k) + 1] = 0.0;
      }
      if (k + 2 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          s = fabs(h[1]) + fabs(h[5]);
          tst = 0.75 * s + h[0];
          aa = -0.4375 * s;
          ab = s;
          h22 = tst;
        } else if (its == 20) {
          s = fabs(h[i + 3 * (i - 1)]) + fabs(h[i - 1]);
          tst = 0.75 * s + h[i + 3 * i];
          aa = -0.4375 * s;
          ab = s;
          h22 = tst;
        } else {
          nr = i + 3 * (i - 1);
          tst = h[nr - 1];
          ab = h[nr];
          aa = h[(i + 3 * i) - 1];
          h22 = h[i + 3 * i];
        }
        s = ((fabs(tst) + fabs(aa)) + fabs(ab)) + fabs(h22);
        if (s == 0.0) {
          rt1r = 0.0;
          tst = 0.0;
          bb = 0.0;
          ab = 0.0;
        } else {
          tst /= s;
          ab /= s;
          aa /= s;
          h22 /= s;
          ba = (tst + h22) / 2.0;
          tst = (tst - ba) * (h22 - ba) - aa * ab;
          ab = sqrt(fabs(tst));
          if (tst >= 0.0) {
            rt1r = ba * s;
            bb = rt1r;
            tst = ab * s;
            ab = -tst;
          } else {
            rt1r = ba + ab;
            bb = ba - ab;
            if (fabs(rt1r - h22) <= fabs(bb - h22)) {
              rt1r *= s;
              bb = rt1r;
            } else {
              bb *= s;
              rt1r = bb;
            }
            tst = 0.0;
            ab = 0.0;
          }
        }
        m = i - 1;
        if (i - 1 >= 1) {
          aa = h[0] - bb;
          s = (fabs(aa) + fabs(ab)) + fabs(h[1]);
          ba = h[1] / s;
          v[0] = (ba * h[3] + (h[0] - rt1r) * (aa / s)) - tst * (ab / s);
          v[1] = ba * (((h[0] + h[4]) - rt1r) - bb);
          v[2] = ba * h[5];
          s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
        }
        for (b_k = m; b_k <= i; b_k++) {
          nr = (i - b_k) + 2;
          if (3 < nr) {
            nr = 3;
          }
          if (b_k > i - 1) {
            hoffset = (b_k + 3 * (b_k - 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }
          tst = v[0];
          ba = xzlarfg(nr, &tst, v);
          v[0] = tst;
          if (b_k > i - 1) {
            h[b_k - 1] = tst;
            h[b_k] = 0.0;
            if (b_k < i) {
              h[2] = 0.0;
            }
          }
          d = v[1];
          ab = ba * v[1];
          if (nr == 3) {
            s = v[2];
            tst = ba * v[2];
            for (j = b_k; j < 4; j++) {
              sum1_tmp_tmp = 3 * (j - 1);
              hoffset = b_k + sum1_tmp_tmp;
              aa = (h[hoffset - 1] + d * h[hoffset]) + s * h[sum1_tmp_tmp + 2];
              h[hoffset - 1] -= aa * ba;
              h[hoffset] -= aa * ab;
              h[sum1_tmp_tmp + 2] -= aa * tst;
            }
            if (b_k + 3 < i + 1) {
              sum1_tmp_tmp = b_k + 2;
            } else {
              sum1_tmp_tmp = i;
            }
            for (j = 0; j <= sum1_tmp_tmp; j++) {
              hoffset = j + 3 * (b_k - 1);
              nr = j + 3 * b_k;
              aa = (h[hoffset] + d * h[nr]) + s * h[j + 6];
              h[hoffset] -= aa * ba;
              h[nr] -= aa * ab;
              h[j + 6] -= aa * tst;
            }
            for (j = 0; j < 3; j++) {
              hoffset = j + 3 * (b_k - 1);
              nr = j + 3 * b_k;
              aa = (z[hoffset] + d * z[nr]) + s * z[j + 6];
              z[hoffset] -= aa * ba;
              z[nr] -= aa * ab;
              z[j + 6] -= aa * tst;
            }
          } else if (nr == 2) {
            for (j = b_k; j < 4; j++) {
              hoffset = b_k + 3 * (j - 1);
              tst = h[hoffset - 1];
              aa = tst + d * h[hoffset];
              h[hoffset - 1] = tst - aa * ba;
              h[hoffset] -= aa * ab;
            }
            for (j = 0; j <= i; j++) {
              hoffset = j + 3 * (b_k - 1);
              nr = j + 3 * b_k;
              aa = h[hoffset] + d * h[nr];
              h[hoffset] -= aa * ba;
              h[nr] -= aa * ab;
            }
            for (j = 0; j < 3; j++) {
              hoffset = j + 3 * (b_k - 1);
              tst = z[hoffset];
              nr = j + 3 * b_k;
              aa = tst + d * z[nr];
              z[hoffset] = tst - aa * ba;
              z[nr] -= aa * ab;
            }
          }
        }
        its++;
      }
    }
    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((k + 2 != i + 1) && (k + 2 == i)) {
        sum1_tmp_tmp = i + 3 * i;
        d = h[sum1_tmp_tmp - 1];
        m = 3 * (i - 1);
        nr = i + m;
        s = h[nr];
        tst = h[sum1_tmp_tmp];
        xdlanv2(&h[(i + 3 * (i - 1)) - 1], &d, &s, &tst, &ab, &aa, &ba, &bb,
                &h22, &rt1r);
        h[sum1_tmp_tmp - 1] = d;
        h[nr] = s;
        h[sum1_tmp_tmp] = tst;
        if (3 > i + 1) {
          nr = 1 - i;
          hoffset = i + (i + 1) * 3;
          for (b_k = 0; b_k <= nr; b_k++) {
            sum1_tmp_tmp = hoffset + b_k * 3;
            tst = h22 * h[sum1_tmp_tmp - 1] + rt1r * h[sum1_tmp_tmp];
            h[sum1_tmp_tmp] =
                h22 * h[sum1_tmp_tmp] - rt1r * h[sum1_tmp_tmp - 1];
            h[sum1_tmp_tmp - 1] = tst;
          }
        }
        nr = i * 3;
        if (i - 1 >= 1) {
          tst = h22 * h[m] + rt1r * h[nr];
          h[nr] = h22 * h[nr] - rt1r * h[m];
          h[m] = tst;
        }
        tst = h22 * z[m] + rt1r * z[nr];
        z[nr] = h22 * z[nr] - rt1r * z[m];
        z[m] = tst;
        tst = z[nr + 1];
        ab = z[m + 1];
        z[nr + 1] = h22 * tst - rt1r * ab;
        z[m + 1] = h22 * ab + rt1r * tst;
        tst = z[nr + 2];
        ab = z[m + 2];
        z[nr + 2] = h22 * tst - rt1r * ab;
        z[m + 2] = h22 * ab + rt1r * tst;
      }
      i = k;
    }
  }
  return info;
}

/* End of code generation (xhseqr.c) */
