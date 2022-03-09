/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlartg.c
 *
 * Code generation for function 'xzlartg'
 *
 */

/* Include files */
#include "xzlartg.h"
#include "PolarDecompositionFast_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn)
{
  double f2;
  double fs_im;
  double fs_re;
  double g2;
  double g2s;
  double gs_im;
  double gs_re;
  double scale;
  double scale_tmp;
  boolean_T guard1 = false;
  scale_tmp = fabs(f.re);
  f2 = fabs(f.im);
  if (f2 > scale_tmp) {
    scale_tmp = f2;
  }
  f2 = fabs(g.re);
  scale = fabs(g.im);
  if (scale > f2) {
    f2 = scale;
  }
  scale = scale_tmp;
  if (f2 > scale_tmp) {
    scale = f2;
  }
  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));
    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));
      guard1 = true;
    }
  } else {
    guard1 = true;
  }
  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0 > g2) {
      scale = 1.0;
    }
    if (f2 <= scale * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        g2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        if (scale_tmp > 1.0) {
          g2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          f2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          g2 = rt_hypotd_snf(f2, scale);
          fs_re = f2 / g2;
          fs_im = scale / g2;
        }
        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      scale = sqrt(g2 / f2 + 1.0);
      *cs = 1.0 / scale;
      g2 += f2;
      fs_re = scale * fs_re / g2;
      fs_im = scale * fs_im / g2;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
             creal_T *r)
{
  double f2;
  double f2s;
  double fs_im;
  double fs_re;
  double g2;
  double gs_im;
  double gs_re;
  double scale;
  double scale_tmp;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  scale_tmp = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale_tmp) {
    scale_tmp = f2s;
  }
  f2s = fabs(g.re);
  scale = fabs(g.im);
  if (scale > f2s) {
    f2s = scale;
  }
  scale = scale_tmp;
  if (f2s > scale_tmp) {
    scale = f2s;
  }
  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));
    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));
      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }
  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0 > g2) {
      scale = 1.0;
    }
    if (f2 <= scale * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        f2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / f2;
        sn->im = -gs_im / f2;
      } else {
        g2 = sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2;
        if (scale_tmp > 1.0) {
          f2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / f2;
          fs_im = f.im / f2;
        } else {
          scale = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          f2 = rt_hypotd_snf(scale, f2s);
          fs_re = scale / f2;
          fs_im = f2s / f2;
        }
        gs_re /= g2;
        gs_im = -gs_im / g2;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = sqrt(g2 / f2 + 1.0);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      f2 += g2;
      f2s = r->re / f2;
      scale = r->im / f2;
      sn->re = f2s * gs_re - scale * -gs_im;
      sn->im = f2s * -gs_im + scale * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else if (rescaledir < 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 1.3435752215134178E-138;
          r->im *= 1.3435752215134178E-138;
        }
      }
    }
  }
}

/* End of code generation (xzlartg.c) */
