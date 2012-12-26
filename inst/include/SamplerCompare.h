// From SamplerCompare, (c) 2010 Madeleine Thompson

#ifndef SAMPLERCOMPARE_H
#define SAMPLERCOMPARE_H

#include <Rinternals.h>

// A structure representing a log-density function.  ndim specifies
// the number of dimensions of the distribution, log_dens specifies the log
// density itself, and context contains ancillary data.  The object
// itself should always be the first parameter passed to log_dens.
// The second is an ndim-long array of doubles indicating a coordinate
// set.  The return value is a double indicating a log density.  If
// the third argument is nonzero, that means the gradient has been
// requested.  In that case, the fourth argument is an ndim-long array
// of doubles, which will be filled in with the gradient of the log
// density at the second argument.

struct dist_t_forward;

typedef double log_density_t(struct dist_t_forward*, double*, int, double*);

typedef struct dist_t_forward {
  log_density_t *log_dens;
  SEXP context;
  int ndim;
} dist_t;

// API for a sampler implemented in C.  x0, the start state, has
// length ds->ndim.  X_out is column-major with dimensions sample_size
// by ds->ndim, and should be filled in by the implementing function.
// The caller must wrap calls to a sampler_t with GetRNGstate() and
// PutRNGstate().

typedef void sampler_t(SEXP sampler_context, dist_t *ds, double *x0,
                       int sample_size, double tuning, double *X_out);

// Alternate API for a sampler implemented in C.  See transition_sample
// in src/slice.c and the section of "glue" vignette entitled "The
// Transition Function API" for more information on how this works.

typedef void transition_fn(SEXP sampler_context, dist_t *ds, double *x0,
                           double tuning, double *x1);

#endif /* SAMPLERCOMPARE_H */
