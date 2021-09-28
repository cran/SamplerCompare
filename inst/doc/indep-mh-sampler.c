// From SamplerCompare, (c) 2010 Madeleine Thompson

// indep-mh-sampler.c demonstrates how to implement a sampler
// in C.  It defines an independence-Metropolis-Hastings sampler.
// See the vignette "R/C Glue in SamplerCompare" for more information
// on the API.  See Markov Chain Monte Carlo in Practice (1995, Gilks,
// Richardson, & Spiegelhalter, eds.), pp. 9-10 for a description of
// the algorithm.

#include <R.h>
#include <R_ext/BLAS.h>

#include <SamplerCompare.h>

void indep_mh(SEXP sampler_context, dist_t *ds, double *x0, int sample_size,
              double tuning, double *X_out) {
  double proposal_mean = REAL(sampler_context)[0];  // extract prop. mean
  double proposal_stddev = tuning;
  double x[ds->ndim];
  int one = 1;
  dcopy_(&ds->ndim, x0, &one, x, &one);     // copy x0 to x
  double y = ds->log_dens(ds, x, 0, NULL);  // log dens. at x

  for (int obs = 0; obs < sample_size; obs++) {
    R_CheckUserInterrupt();  // so user can ^C
    for (int dim = 0; dim < ds->ndim; dim++) {
      double x_dim_old = x[dim];                               // stash x[dim]
      x[dim] = proposal_mean + proposal_stddev * norm_rand();  // propose new x
      double y_prop = ds->log_dens(ds, x, 0, NULL);  // log dens. at prop.
      if (unif_rand() > exp(y_prop - y))             // If we reject prop.,
        x[dim] = x_dim_old;                          // restore x[dim].
      else                                           // Otherwise, store
        y = y_prop;                                  // log dens for new x.
    }
    dcopy_(&ds->ndim, x, &one, X_out + obs, &sample_size);  // record state
  }
}
