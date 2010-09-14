// From SamplerCompare, (c) 2010 Madeleine Thompson
// $Id: slice.c 1521 2010-08-28 21:17:53Z mthompson $

// slice.c contains a C implementation of the shrinking rank slice
// sampler.  It uses the transition_fn abstraction that allows a sampler
// to be constructed from a group of functions that leave the target
// distribution invariant.  This same abstraction is used by several
// other samplers that are not a part of the SamplerCompare package.

#include <R.h>
#include <R_ext/BLAS.h>

#include "../inst/include/SamplerCompare.h"

static void projv(double *J, int ncol, int nrow, double *v);

sampler_t transition_sample;
transition_fn sr_draw;

// transition_sample is a generic MCMC sampler wrapper parameterized
// on a function that transitions from state to state.  Its context
// is either an R raw with a function pointer with a signature of type
// transition_fn, or a list where the first element is such a raw.
// If it is a list, the remaining elements of the list can be used
// for tuning parameters for the transition function itself.  See
// "The Transition Function API" in the "glue" vignette for more
// information.

void transition_sample(SEXP context, dist_t *ds, double *x0,
                       int sample_size, double tuning, double *X_out) {

  // Extract the transition function.

  transition_fn *T = NULL;
  if (TYPEOF(context)==RAWSXP)
    T = *(transition_fn**)RAW(context);
  else if (isNewList(context) && TYPEOF(VECTOR_ELT(context,0))==RAWSXP)
    T = *(transition_fn**)RAW(VECTOR_ELT(context, 0));
  else
    error("transition_sample context must be either a raw function pointer "
          "or a list with a raw function pointer as its first element\n");

  // Copy the initial state to x and X_out.

  int one = 1;
  double x[ds->ndim];
  dcopy_(&ds->ndim, x0, &one, x, &one);
  dcopy_(&ds->ndim, x0, &one, X_out, &sample_size);

  // Call the transition function sample_size times, storing the
  // resulting state in X_out after each call.

  for (int obs=1; obs<sample_size; obs++) {
    R_CheckUserInterrupt();
    double x1[ds->ndim];
    T(context, ds, x, tuning, x1);
    dcopy_(&ds->ndim, x1, &one, x, &one);
    dcopy_(&ds->ndim, x1, &one, X_out+obs, &sample_size);
  }
}

// A transition function for implementing shrinking rank slice
// sampling (see "Covariance Adaptive Slice Sampling," Thompson and
// Neal, 2010, sec. 5).  sigma_c, the "tuning" parameter, represents
// the initial crumb standard deviation.  The scaling factor, which
// is usually left constant at 0.95, can be passed as the second element
// of the sampler context, where the first is the function pointer to
// sr_draw itself, used by transition_sample().
//
// The context may be set to an R NULL if the caller wants to use
// this as part of a compound sampling operation and does not have
// a convenient sampler context in the right form.  In that case, the
// scaling factor is fixed at 0.95.
//
// See the R help (?shrinking.rank.sample) for more information on
// the algorithm itself.

void sr_draw(SEXP sampler_context, dist_t *ds, double *x0,
             double sigma_c, double *x1) {

  // Identify a scaling factor and minimum dimension.

  double downscale = 0.95;
  int min_dim = 1;
  if (!isNull(sampler_context)) {
    downscale = *REAL(VECTOR_ELT(sampler_context, 1));
    min_dim = *INTEGER(VECTOR_ELT(sampler_context, 2));
  }
  if (!(downscale>0) || !(downscale<=1))
    error("invalid scaling factor: %.3g", downscale);

  // Choose a slice level.

  double y_slice = ds->log_dens(ds, x0, 0, NULL) - exp_rand(); // slice level

  // J needs to be allocated dynamically since many calls to
  // sr_draw may occur between times control returns to R and
  // R_alloc-allocated memory is freed.

  double *J = NULL;
  int ncol_J = 0;
  int ncrumb = 0;

  // Set cbar, the precision-weighted sum of crumbs so far
  // drawn, to zero.  prec_cbar, the associated precision, will
  // reflect the infinite posterior uncertainty with an initial
  // value of zero.  cbar is maintained as an offset to x0, not
  // the origin.

  double cbar[ds->ndim], d_zero = 0.0, d_one = 1.0;
  const int zero = 0, one = 1;
  dcopy_(&ds->ndim, &d_zero, &zero, cbar, &one);
  double prec_cbar = 0.0;

  // Initialize the crumb standard deviation to sigma_c.

  double sigma_ck = sigma_c;

  while(1) {
    // Draw a new crumb, c, with standard deviation sigma_ck.  Like
    // cbar, the crumb is stored relative to x0.

    ncrumb++;
    double c[ds->ndim];
    for (int i=0; i<ds->ndim; i++)
      c[i] = norm_rand();
    dscal_(&ds->ndim, &sigma_ck, c, &one);

    // Update cbar and prec_cbar to reflect the newly drawn crumb.
    // cbar takes the value cbar+c/sigma_ck^2, and prec_cbar increases
    // sigma_ck^(-2).

    double prec_c = 1.0/(sigma_ck*sigma_ck);
    daxpy_(&ds->ndim, &prec_c, c, &one, cbar, &one);
    prec_cbar = prec_cbar + prec_c;

    // Draw a proposal as:
    //
    // x1 = x0 + projv(J, prec_cbar^(-1) * cbar + prec_cbar^(-1/2) * rnorm(p))
    //
    // x0 + projv(J, cbar/prec_cbar) is the posterior mean of x0, and
    // 1/sqrt(prec_cbar) is its posterior variance, given the crumbs.

    for (int i=0; i<ds->ndim; i++)
      x1[i] = norm_rand();
    double scale = 1.0/sqrt(prec_cbar);
    dscal_(&ds->ndim, &scale, x1, &one);
    scale = 1.0/prec_cbar;
    daxpy_(&ds->ndim, &scale, cbar, &one, x1, &one);
    projv(J, ncol_J, ds->ndim, x1);
    daxpy_(&ds->ndim, &d_one, x0, &one, x1, &one);

    // Check to see if the proposal is in the slice.  If so accept it.

    double grad[ds->ndim], y1;
    y1 = ds->log_dens(ds, x1, 1, grad);
    if (y1 >= y_slice)
      break;

    // If we are outside the support of the distribution, our
    // proposal distribution is probably much too large.  Shrink by
    // an order of magnitude and try again.

    if (!R_FINITE(y1)) {
      sigma_ck = 0.1 * downscale * sigma_ck;
      continue;
    }

    // If a column is not added to J because it is already high-rank,
    // downscale the crumb standard deviation and try again.

    if (ds->ndim - ncol_J <= min_dim) {
      sigma_ck = downscale * sigma_ck;
      continue;
    }

    // Now, we consider adding a column to J.  Start by projecting
    // the gradient into the nullspace of J and normalizing it.

    projv(J, ncol_J, ds->ndim, grad);
    double inv_normg = 1.0/dnrm2_(&ds->ndim, grad, &one);
    dscal_(&ds->ndim, &inv_normg, grad, &one);

    // Compute angle of projected gradient with the gradient.  If
    // it is less than 60 deg., add the projected gradient as a
    // column of J.  This should handle infinities gracefully (by
    // failing to add the column).  If a column is not added to
    // J because the angle is too large, downscale the crumb
    // standard deviation.

    double cos_theta = ddot_(&ds->ndim, grad, &one, grad, &one) /
      dnrm2_(&ds->ndim, grad, &one);
    if (cos_theta > 0.5) {
      J = Realloc(J, ds->ndim*(ncol_J+1), double);
      dcopy_(&ds->ndim, grad, &one, J+(ds->ndim)*ncol_J, &one);
      ncol_J++;
    } else {
      sigma_ck = downscale * sigma_ck;
    }
  }

  if (J!=NULL)
    Free(J);
}

// Project the vector v into the nullspace of the matrix J.  J is
// an nrow*ncol matrix stored in column major order; v is an nrow-length
// vector.
//
// In R notation, the operation performed is:
//
//   v <- J %*% t(J) %*% v
//
// If J has no columns, this function does nothing.

static void projv(double *J, int ncol, int nrow, double *v) {
  if (ncol==0)
    return;
  double v2[ncol], d_zero = 0.0, d_one = 1.0, d_minus_one = -1.0;
  const int one = 1, zero = 0;
  dcopy_(&ncol, &d_zero, &zero, v2, &one);
  dgemv_("t", &nrow, &ncol, &d_one, J, &nrow, v, &one, &d_zero, v2, &one);
  dgemv_("n", &nrow, &ncol, &d_minus_one, J, &nrow, v2, &one, &d_one, v, &one);
}
