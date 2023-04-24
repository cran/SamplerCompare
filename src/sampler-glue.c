// From SamplerCompare, (c) 2010 Madeleine Thompson

// sampler-glue.c contains functions related to the R/C interface
// for samplers and distributions.  See the vignette "R/C Glue in
// SamplerCompare" for more information.

#define USE_FC_LEN_T

#include <stdint.h>
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include "../inst/include/SamplerCompare.h"
#include "init.h"

// structure representing a function passed into a sampler, to be used
// in a log_density_t.  log_dens_fun is the R function object, envir
// is the environment to evaluate it in, and calls is the number of
// times the function has been called from within the sampler.

typedef struct {
  SEXP log_dens;
  SEXP envir;
  int evals;
  int grads;
} R_stub_context_t;

// similar to R_stub_context_t, just with information for calling
// a C function instead of an R function.

typedef struct {
  dist_t ds;  // the inner dist_t to be invoked
  int evals;
  int grads;
} C_stub_context_t;

static log_density_t C_log_density_stub_func, R_log_density_stub_func;
static SEXP void_as_raw(void *p);
static inline void *raw_as_void(SEXP raw_vec) { return *(void **)RAW(raw_vec); }

// Use a sampler implemented in C to draw a sample from a distribution
// defined in R.  The arguments are:
//
//   sampler          name of symbol of the log density/gradient function
//   sampler_context  an opaque context passed to sampler
//   log_dens         a SEXP of a log.density.and.grad function
//                    implemented in R (see make.dist for convention)
//   x0               a numeric vector containing an initial point
//   sample_size      number of observations requested
//   tuning           one-element numeric vector specifying a tuning parameter
//   envir            an environment to evaluate log_dens in
//
// The return value is an R list with elements:
//
//   X      a matrix containing the rows of observations
//   evals  an integer indicating the number of times log_dens was called
//   grads  an integer indicating the number of times log_dens was
//          called with a gradient requested

SEXP sampler_glue_R_dist(SEXP sampler, SEXP sampler_context, SEXP log_dens,
                         SEXP x0, SEXP sample_size, SEXP tuning, SEXP envir) {
  // Check parameters for validity and unpack some of them into C types.

  if (!isEnvironment(envir)) {
    error("envir is not an environment.");
  }

  int sample_size_int = asInteger(sample_size);
  if (sample_size_int < 1) {
    error("sample size must be a positive integer.");
  }
  int ndim = length(x0);

  double tuning_dbl = asReal(tuning);
  double *x0_dbl = REAL(x0);

  // Locate the sampler as a function pointer.

  if (!isString(sampler)) {
    error("sampler is not a character string.");
  }
  sampler_t *sampler_fp =
      (sampler_t *)R_FindSymbol(CHAR(STRING_ELT(sampler, 0)), "", NULL);
  if (sampler_fp == NULL) {
    error("Cannot locate symbol \"%s\".", CHAR(STRING_ELT(sampler, 0)));
  }

  // Create a stub for log_dens so that it looks like a C density
  // to the sampler.

  R_stub_context_t stub_context = {
      .log_dens = log_dens, .envir = envir, .evals = 0, .grads = 0};
  SEXP raw_context;
  PROTECT(raw_context = void_as_raw(&stub_context));
  dist_t stub_ds = {.log_dens = R_log_density_stub_func,
                    .context = raw_context,
                    .ndim = ndim};

  // Allocate a result matrix, set up the RNG, and call the sampler
  // to draw a sample.

  SEXP X_out;
  PROTECT(X_out = allocMatrix(REALSXP, sample_size_int, ndim));
  GetRNGstate();
  sampler_fp(sampler_context, &stub_ds, x0_dbl, sample_size_int, tuning_dbl,
             REAL(X_out));
  PutRNGstate();

  // Set up return value as an R object.

  SEXP ans, ans_names;
  PROTECT(ans = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ans, 0, X_out);
  SET_VECTOR_ELT(ans, 1, ScalarInteger(stub_context.evals));
  SET_VECTOR_ELT(ans, 2, ScalarInteger(stub_context.grads));
  PROTECT(ans_names = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ans_names, 0, mkString("X"));
  SET_VECTOR_ELT(ans_names, 1, mkString("evals"));
  SET_VECTOR_ELT(ans_names, 2, mkString("grads"));
  setAttrib(ans, R_NamesSymbol, ans_names);
  UNPROTECT(4);

  return ans;
}

// This function wraps an R log density function so that it exposes
// the interface expected by a sampler_t and keeps track of the number
// of times it is called.

static double R_log_density_stub_func(dist_t *ds, double *x, int compute_grad,
                                      double *grad) {
  SEXP xsexp, fcall, result_sexp, compute_grad_sexp, result_names;

  R_stub_context_t *stub_context = (R_stub_context_t *)raw_as_void(ds->context);

  // Allocate R variables for the arguments to the R log.density.and.grad
  // and call it.

  PROTECT(xsexp = allocVector(REALSXP, ds->ndim));
  memmove(REAL(xsexp), x, sizeof(double) * ds->ndim);
  PROTECT(compute_grad_sexp = allocVector(LGLSXP, 1));
  LOGICAL(compute_grad_sexp)[0] = (compute_grad != 0);
  PROTECT(fcall = lang3(stub_context->log_dens, xsexp, compute_grad_sexp));
  PROTECT(result_sexp = eval(fcall, stub_context->envir));

  double log_dens = NAN;
  int found_log_dens = 0, found_grad = 0;

  // Unpack the results from the log.density.and.grad into the
  // variable log_dens and (if appropriate) the memory pointed to by
  // grad.

  if (!isNewList(result_sexp)) {
    error("log density function must return a list.");
  }
  PROTECT(result_names = getAttrib(result_sexp, R_NamesSymbol));
  for (int i = 0; i < length(result_sexp); i++) {
    if (!strcmp(CHAR(STRING_ELT(result_names, i)), "log.density")) {
      log_dens = asReal(VECTOR_ELT(result_sexp, i));
      found_log_dens = 1;
    }
    if (compute_grad &&
        !strcmp(CHAR(STRING_ELT(result_names, i)), "grad.log.density")) {
      memmove(grad, REAL(VECTOR_ELT(result_sexp, i)),
              sizeof(double) * ds->ndim);
      found_grad = 1;
    }
  }

  UNPROTECT(5);

  // Throw an error if the log density did not return the appropriate
  // list elements.

  if (!found_log_dens) {
    error("log density did not return log.density element.");
  }
  if (!found_grad && compute_grad) {
    error("log density did not return grad.log.density element.");
  }

  // Increment the evaluation counters.

  stub_context->evals++;
  if (compute_grad) {
    stub_context->grads++;
  }

  return log_dens;
}

// Use a sampler implemented in C to draw a sample from a distribution
// implemented in C.  It would be possible to call the distribution
// through its R glue, but having this function allows C samplers to
// avoid a round trip through the R interpreter.  This is called from
// a wrapper defined in wrap.c.sampler.  The arguments are:
//
//   sampler_name     name of symbol of a sampler_t defining the sampler.
//   sampler_context  an opaque context passed to the sampler
//   log_dens_name    name of symbol of a log_density_t defining
//                    the distribution.
//   dist_context     opaque context passwd to the log density.
//   x0               a numeric vector containing an initial point
//   sample_size      number of observations requested
//   tuning           one-element numeric vector specifying a tuning parameter
//
// The return value is an R list with elements:
//
//   X      a matrix containing the rows of observations
//   evals  an integer indicating the number of times log density was called
//   grads  an integer indicating the number of times log density was
//          called with a gradient requested

SEXP sampler_glue_C_dist(SEXP sampler_name, SEXP sampler_context,
                         SEXP log_dens_name, SEXP dist_context, SEXP x0,
                         SEXP sample_size, SEXP tuning) {
  // Locate symbol for sampler function.

  const char *sampler_str = CHAR(STRING_ELT(sampler_name, 0));
  sampler_t *sampler_fp = (sampler_t *)R_FindSymbol(sampler_str, "", NULL);
  if (sampler_fp == NULL) {
    error("Cannot locate symbol \"%s\".", sampler_str);
  }

  // Locate symbol for log density.

  const char *log_dens_str = CHAR(STRING_ELT(log_dens_name, 0));
  log_density_t *log_dens_fp =
      (log_density_t *)R_FindSymbol(log_dens_str, "", NULL);
  if (log_dens_fp == NULL) {
    error("Cannot locate symbol \"%s\".", log_dens_str);
  }

  // Define a stub function to keep track of the number of function calls.

  int ndim = length(x0);
  C_stub_context_t stub_context = {
      .ds = {.log_dens = log_dens_fp, .ndim = ndim, .context = dist_context},
      .evals = 0,
      .grads = 0};
  SEXP raw_context;
  PROTECT(raw_context = void_as_raw(&stub_context));
  dist_t stub_ds = {.log_dens = C_log_density_stub_func,
                    .context = raw_context,
                    .ndim = ndim};

  // Create a matrix to store the states in and call the sampler.

  SEXP X;
  PROTECT(X = allocMatrix(REALSXP, *REAL(sample_size), ndim));
  GetRNGstate();
  sampler_fp(sampler_context, &stub_ds, REAL(x0), *REAL(sample_size),
             *REAL(tuning), REAL(X));
  PutRNGstate();

  // Construct the result to return.

  const char *result_names[] = {"X", "evals", "grads", ""};
  SEXP result;
  PROTECT(result = mkNamed(VECSXP, result_names));
  SET_VECTOR_ELT(result, 0, X);
  SET_VECTOR_ELT(result, 1, ScalarInteger(stub_context.evals));
  SET_VECTOR_ELT(result, 2, ScalarInteger(stub_context.grads));
  UNPROTECT(3);
  return result;
}

// A wrapper around a log density used by sampler_glue_C_dist to
// count the number of times the log density is called.  Obeys the
// calling convention of log_density_t.

static double C_log_density_stub_func(dist_t *ds, double *x, int compute_grad,
                                      double *grad) {
  C_stub_context_t *stub_context = (C_stub_context_t *)raw_as_void(ds->context);
  stub_context->evals++;
  if (compute_grad) stub_context->grads++;
  return stub_context->ds.log_dens(&stub_context->ds, x, compute_grad, grad);
}

// See ?raw_symbol and the 'glue' vignette for more information.

SEXP raw_symbol(SEXP symbol_name) {
  // Find a function pointer for the requested symbol.

  if (!isString(symbol_name) || length(symbol_name) != 1) {
    error("Invalid symbol_name.");
  }
  const char *symbol_char = CHAR(STRING_ELT(symbol_name, 0));
  void *symbol = (void*)(intptr_t)R_FindSymbol(symbol_char, "", NULL);
  if (symbol == NULL) {
    error("Could not locate symbol \"%s\".", symbol_char);
  }

  // Copy the function pointer to a raw vector and return it.

  return void_as_raw(symbol);
}

// Allocates a raw vector and stores a void* in it.  This allows C
// code to pass around data inside R objects.  The inverse operation
// is raw_as_void.

static SEXP void_as_raw(void *p) {
  SEXP raw_vec;
  PROTECT(raw_vec = allocVector(RAWSXP, sizeof(void *)));
  *(void **)RAW(raw_vec) = p;
  UNPROTECT(1);
  return raw_vec;
}

// This function wraps a log density implemented in C so that it can
// be called from R.  It is used by make.c.dist.  The arguments are:
//
//   c_sym          an R raw vector containing a function pointer
//                  obeying the sampler_t convention, representing
//                  the target distribution's log density.  Often
//                  generated by raw_symbol().
//
//   context        an opaque context passed to the log density
//
//   x              a numeric vector indicating the point at which
//                  to evaluate the log density
//
//   compute_grad   a one-element logical specifying whether the
//                  gradient is requested
//
// The return value is list containing elements "log.density" and
// optionally "grad.log.density".  log.density is a numeric scalar.
// If compute_grad is false, grad.log.density is not set.  Otherwise,
// it is a numeric vector with the same number of elements as x.

SEXP R_invoked_C_glue(SEXP c_sym, SEXP context, SEXP x, SEXP compute_grad) {
  // Identify the problem dimension and extract a function pointer
  // for the log density.

  int ndim = length(x);
  log_density_t *log_dens_fn = *(log_density_t **)RAW(c_sym);
  dist_t ds = {.log_dens = log_dens_fn, .context = context, .ndim = ndim};

  if (*LOGICAL(compute_grad)) {
    // Allocate a vector for the gradient and call the log density.

    SEXP grad, result;
    PROTECT(grad = allocVector(REALSXP, ndim));
    double y = log_dens_fn(&ds, REAL(x), 1, REAL(grad));

    // Store the log density and gradient in a list and return it.

    const char *result_names[] = {"log.density", "grad.log.density", ""};
    PROTECT(result = mkNamed(VECSXP, result_names));
    SET_VECTOR_ELT(result, 0, ScalarReal(y));
    SET_VECTOR_ELT(result, 1, grad);

    UNPROTECT(2);
    return result;

  } else {
    // Call the log density function.

    double y = log_dens_fn(&ds, REAL(x), 0, NULL);

    // Store the log density in a list and return it.

    SEXP result;
    const char *result_names[] = {"log.density", ""};
    PROTECT(result = mkNamed(VECSXP, result_names));
    SET_VECTOR_ELT(result, 0, ScalarReal(y));
    UNPROTECT(1);
    return result;
  }
}
