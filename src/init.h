/* Forward-declare internal functions for R_init_SamplerCompare. */

#ifndef SAMPLERCOMPARE_INIT_H
#define SAMPLERCOMPARE_INIT_H

#include <Rinternals.h>

#include "../inst/include/SamplerCompare.h"

SEXP raw_symbol(SEXP);
SEXP R_invoked_C_glue(SEXP c_sym, SEXP context, SEXP x, SEXP compute_grad);
SEXP sampler_glue_R_dist(SEXP sampler, SEXP sampler_context, SEXP log_dens,
                         SEXP x0, SEXP sample_size, SEXP tuning, SEXP envir);
SEXP sampler_glue_C_dist(SEXP sampler_name, SEXP sampler_context,
                         SEXP log_dens_name, SEXP dist_context, SEXP x0,
                         SEXP sample_size, SEXP tuning);

void sr_draw(SEXP sampler_context, dist_t *ds, double *x0,
             double sigma_c, double *x1);
void arms_sample(SEXP sampler_context, dist_t *ds, double *x0, int sample_size,
                 double tuning, double *X_out);
void transition_sample(SEXP context, dist_t *ds, double *x0, int sample_size,
                       double tuning, double *X_out);

double Gauss2_log_dens(dist_t *dist, double *x, int compute_grad, double *grad);
double cone_log_dens(dist_t *dist, double *x, int compute_grad, double *grad);

extern void dchud_(void);
extern void dchdd_(void);

#endif /* SAMPLERCOMPARE_INIT_H */
