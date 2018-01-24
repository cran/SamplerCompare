// From SamplerCompare, (c) 2010 Madeleine Thompson

// distributions.c contains C implementations of a two dimensional
// Gaussian and Roberts and Rosenthal's cone distribution.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <math.h>

#include "../inst/include/SamplerCompare.h"
#include "init.h"

log_density_t Gauss2_log_dens, cone_log_dens;

// This computes a log density and gradient of a two dimensional
// Gaussian.  The context parameter is a pointer to an array of three
// doubles, the first two of which are the mean and the third of which
// is the correlation.  This is intended to test the R-C glue for
// distributions.  The following two commands produce similar
// distribution objects:
//
//   wrap.c.dist(2, 'Gauss2-C', 'Gauss2_fdf', c(1, 2, 0.8), mean=c(1,2))
//   make.gaussian(c(1,2), rho=0.8)
//
// It is only used in test code and is not exported to users.

double Gauss2_log_dens(dist_t *dist, double *x,
                       int compute_grad, double *grad) {
  if (dist->ndim!=2)
    error("Gauss2_log_dens: dimension other than two (%d)", dist->ndim);

  double mu1 = REAL(dist->context)[0];
  double mu2 = REAL(dist->context)[1];
  double rho = REAL(dist->context)[2];
  if (!(rho<=1) || !(rho>=-1))
    error("Gauss2_fdf: invalid correlation (%.6g)", rho);
  double det = 1-rho*rho;

  double z0 = x[0] - mu1;
  double z1 = x[1] - mu2;

  double y = -1.83787706640935 - 0.5 * log(det)
    - 0.5 * (z0*z0 - 2 * rho * z0 * z1 + z1*z1) / det;
  if (compute_grad) {
    grad[0] = -1.0/det * z0 + rho / det * z1;
    grad[1] = -1.0/det * z1 + rho / det * z0;
  }
  return y;
}

// Log density for Roberts and Rosenthal's cone distribution.  See
// ?make.cone.dist for more information.

double cone_log_dens(dist_t *dist, double *x,
                     int compute_grad, double *grad) {
  int one = 1;
  double y = -dnrm2_(&dist->ndim, x, &one);
  if (compute_grad) {
    for (int i=0; i<dist->ndim; i++)
      grad[i] = -fabs(x[i]);
  }
  return(y);
}
