#include "bessel.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>

#include "const.h"

typedef struct {
  BfReal * c;   /* coefficients                */
  int order;    /* order of expansion          */
  BfReal a;     /* lower interval point        */
  BfReal b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
} cheb_series;

static BfReal cheb_eval(const cheb_series * cs, BfReal x) {
  BfReal d  = 0.0;
  BfReal dd = 0.0;

  BfReal y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  BfReal y2 = 2.0 * y;

  for (int j = cs->order; j >= 1; j--) {
    BfReal temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  d = y*d - dd + 0.5 * cs->c[0];

  return d;
}

static BfReal bj0_data[13] = {
   0.100254161968939137,
  -0.665223007764405132,
   0.248983703498281314,
  -0.0332527231700357697,
   0.0023114179304694015,
  -0.0000991127741995080,
   0.0000028916708643998,
  -0.0000000612108586630,
   0.0000000009838650793,
  -0.0000000000124235515,
   0.0000000000001265433,
  -0.0000000000000010619,
   0.0000000000000000074,
};

static cheb_series bj0_cs = {bj0_data, 12, -1, 1, 9};

static BfReal by0_data[13] = {
  -0.011277839392865573,
  -0.128345237560420350,
  -0.104378847997942490,
   0.023662749183969695,
  -0.002090391647700486,
   0.000103975453939057,
  -0.000003369747162423,
   0.000000077293842676,
  -0.000000001324976772,
   0.000000000017648232,
  -0.000000000000188105,
   0.000000000000001641,
  -0.000000000000000011
};

static cheb_series by0_cs = {by0_data, 12, -1, 1, 8};

static BfReal bm0_data[21] = {
   0.09284961637381644,
  -0.00142987707403484,
   0.00002830579271257,
  -0.00000143300611424,
   0.00000012028628046,
  -0.00000001397113013,
   0.00000000204076188,
  -0.00000000035399669,
   0.00000000007024759,
  -0.00000000001554107,
   0.00000000000376226,
  -0.00000000000098282,
   0.00000000000027408,
  -0.00000000000008091,
   0.00000000000002511,
  -0.00000000000000814,
   0.00000000000000275,
  -0.00000000000000096,
   0.00000000000000034,
  -0.00000000000000012,
   0.00000000000000004
};

static const cheb_series _amp_phase_bm0_cs = {bm0_data, 20, -1, 1, 10};

static BfReal bth0_data[24] = {
  -0.24639163774300119,
   0.001737098307508963,
  -0.000062183633402968,
   0.000004368050165742,
  -0.000000456093019869,
   0.000000062197400101,
  -0.000000010300442889,
   0.000000001979526776,
  -0.000000000428198396,
   0.000000000102035840,
  -0.000000000026363898,
   0.000000000007297935,
  -0.000000000002144188,
   0.000000000000663693,
  -0.000000000000215126,
   0.000000000000072659,
  -0.000000000000025465,
   0.000000000000009229,
  -0.000000000000003448,
   0.000000000000001325,
  -0.000000000000000522,
   0.000000000000000210,
  -0.000000000000000087,
   0.000000000000000036
};

static const cheb_series _amp_phase_bth0_cs = {bth0_data, 23, -1, 1, 12};

static BfReal cos_pi4_plus_eps(BfReal y, BfReal eps) {
  BfReal sy = sin(y);
  BfReal cy = cos(y);
  BfReal s = sy + cy;
  BfReal d = sy - cy;
  BfReal seps;
  BfReal ceps;
  if(fabs(eps) < BF_ROOT5_EPS) {
    BfReal e2 = eps*eps;
    seps = eps * (1.0 - e2/6.0 * (1.0 - e2/20.0));
    ceps = 1.0 - e2/2.0 * (1.0 - e2/12.0);
  }
  else {
    seps = sin(eps);
    ceps = cos(eps);
  }
  return (ceps * s - seps * d)/ BF_SQRT2;
}

BfReal bf_j0(BfReal x) {
  BfReal y = fabs(x);

  if (y < 2.0*BF_SQRT_EPS)
    return 1.0;

  if (y <= 4.0)
    return cheb_eval(&bj0_cs, 0.125*y*y - 1.0);

  BfReal z = 32.0/(y*y) - 1.0;
  BfReal ca_val = cheb_eval(&_amp_phase_bm0_cs, z);
  BfReal ct_val = cheb_eval(&_amp_phase_bth0_cs, z);
  BfReal cp_val = cos_pi4_plus_eps(y, ct_val/y);
  BfReal sqrty = sqrt(y);
  BfReal ampl  = (0.75 + ca_val) / sqrty;
  return ampl * cp_val;
}

static BfReal sin_pi4_plus_eps(BfReal y, BfReal eps) {
  BfReal sy = sin(y);
  BfReal cy = cos(y);
  BfReal s = sy + cy;
  BfReal d = sy - cy;
  BfReal seps;
  BfReal ceps;
  if(fabs(eps) < BF_ROOT5_EPS) {
    BfReal e2 = eps*eps;
    seps = eps * (1.0 - e2/6.0 * (1.0 - e2/20.0));
    ceps = 1.0 - e2/2.0 * (1.0 - e2/12.0);
  } else {
    seps = sin(eps);
    ceps = cos(eps);
  }
  return (ceps * d + seps * s)/ BF_SQRT2;
}

BfReal bf_y0(BfReal x) {
  BfReal two_over_pi = 2.0/BF_PI;
  BfReal xmax        = 1.0/BF_EPS;

  assert(x > 0);

  if(x < 4.0) {
    BfReal J0 = bf_j0(x);
    BfReal c_val = cheb_eval(&by0_cs, 0.125*x*x-1.0);
    return two_over_pi*(-BF_LN2 + log(x))*J0 + 0.375 + c_val;
  }

  if(x < xmax) {
    /* Leading behaviour of phase is x, which is exact,
     * so the error is bounded.
     */
    BfReal z  = 32.0/(x*x) - 1.0;
    BfReal c1_val = cheb_eval(&_amp_phase_bm0_cs, z);
    BfReal c2_val = cheb_eval(&_amp_phase_bth0_cs, z);
    BfReal sp_val = sin_pi4_plus_eps(x, c2_val/x);
    BfReal sqrtx = sqrt(x);
    BfReal ampl  = (0.75 + c1_val) / sqrtx;
    return ampl * sp_val;
  }

  // UNDERFLOW_ERROR(result);
  assert(false); // TODO: throw error
}