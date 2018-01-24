// From SamplerCompare, (c) 2010 Madeleine Thompson

// arms.c contains an implementation of Gilks et al's Adaptive
// Rejection Metropolis Sampling.

// This file uses R's Calloc/Realloc/Free allocator for everything
// except the data to be returned, but does not free anything on
// interrupt.  Since a relatively small amount of working memory is
// allocated, this is probably okay.  Fixing this would introduce more
// complexity than it's worth.

#include <string.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>

#include "../inst/include/SamplerCompare.h"
#include "init.h"

// A structure representing an approximate (exact when unimodal)
// bound on a one dimensional cut through a log density.  There are
// nabsc points, with abscissae in absc and ordinates in logdens, each
// of which is an array nabsc long.  This can be freed with free_envelope.

typedef struct {
  int nabsc;
  double *logdens;
  double *absc;
} envelope_t;

// A structure representing a line segment.  A set of these is
// generated from an envelope, and used to represent a piecewise
// exponential distribution.

typedef struct { double left, right, slope, icept, weight; } lineseg_t;

// These functions are documented before their definitions.

sampler_t arms_sample;
static void free_envelope(envelope_t *envelope);
static void add_to_envelope(envelope_t *envelope, double new_absc,
                            double new_dens);
static double sample_piecewise_exponential(int nlines, lineseg_t *linesegs,
                                           int *degen);
static double sample_chopped_exponential(lineseg_t *ls);
static void envelope_segments(envelope_t *envelope, lineseg_t **linesegs_ptr,
                              int *nlines_ptr);

static envelope_t *identify_bounds(dist_t *ds, double *x, int dim, double w);

static inline double min(double a, double b) { return (a < b) ? a : b; }
static inline double max(double a, double b) { return (a > b) ? a : b; }

static void arms_draw(dist_t *ds, double *x0, double *x1,
                      envelope_t **envelopes, int *rejections);
static double arms_sample1(dist_t *ds, double *x0, int which_dim,
                           envelope_t *envelope, int *rejections);
static void envelope_line(envelope_t *envelope, int i, double *m, double *b);
static double line_intersection(double m1, double b1, double m2, double b2);
static double eval_segment(int nlines, lineseg_t *linesegs, double x);
static void expand_envelope(envelope_t *env, dist_t *ds, double *x0,
                            int which_dim);
static envelope_t *copy_envelope(envelope_t *env);

static void find_decrease(dist_t *ds, double *x, int dim, double w, double y,
                          double *xdec, double *ydec);

static void bsearch_function(dist_t *ds, double *x, int dim, double lower,
                             double upper, double ylower, double yupper,
                             double *q, double *y, int upper_bound);

static const int one = 1;  // for BLAS

static inline void Rassert(int v) {
  if (!v) error("Assertion failed in arms.c.");
}

void arms_sample(SEXP sampler_context, dist_t *ds, double *x0, int sample_size,
                 double tuning, double *X_out) {
  // Initialize envelopes

  envelope_t **envelopes = Calloc(ds->ndim, envelope_t *);
  for (int dim = 0; dim < ds->ndim; dim++)
    envelopes[dim] = identify_bounds(ds, x0, dim, tuning);

  // Use arms_draw sample_size times to draw observations.

  double x[ds->ndim];
  memmove(x, x0, sizeof(double) * ds->ndim);
  int rejections = 0;

  for (int obs = 0; obs < sample_size; obs++) {
    arms_draw(ds, x, x, envelopes, &rejections);
    dcopy_(&ds->ndim, x, &one, &X_out[obs], &sample_size);
  }

  for (int dim = 0; dim < ds->ndim; dim++) {
    free_envelope(envelopes[dim]);
  }
  Free(envelopes);
}

// Uses ARMS to draw a single observation from a distribution, ds,
// starting from state x0.  Stores resulting state in x1.  envelopes
// is an array of initial envelopes; it is not modified.  It may not
// depend on x0.  rejections is a pointer to a counter of M-H rejections.
// If not null, it is incremented every time a proposal is rejected.

static void arms_draw(dist_t *ds, double *x0, double *x1,
                      envelope_t **envelopes, int *rejections) {
  dcopy_(&ds->ndim, x0, &one, x1, &one);

  // Update each coordinate of x0 in turn...

  for (int i = 0; i < ds->ndim; i++) {
    envelope_t *E = copy_envelope(envelopes[i]);
    x1[i] = arms_sample1(ds, x1, i, E, rejections);
    free_envelope(E);
  }
}

// Update a single coordinate in x0 with ARMS, returning the new value.
// ds, x0, and rejections are as with arms_draw.  which_dim
// specifies the coordinate to update.  envelope is an initial guess
// at an envelope, but does not need to satisfy any constraints other
// than finiteness.

static double arms_sample1(dist_t *ds, double *x0, int which_dim,
                           envelope_t *envelope, int *rejections) {
  double x1[ds->ndim];
  memmove(x1, x0, ds->ndim * sizeof(double));
  double y1;
  int nlines;
  lineseg_t *linesegs = NULL;
  while (1) {
    R_CheckUserInterrupt();

    if (linesegs != NULL) {
      Free(linesegs);
    }

    // Ensure the envelope satisfies the slope constraints, transform
    // it into a set of segments, and draw from the implied distribution.

    expand_envelope(envelope, ds, x1, which_dim);
    envelope_segments(envelope, &linesegs, &nlines);
    int degen;
    x1[which_dim] = sample_piecewise_exponential(nlines, linesegs, &degen);

    // Accept the proposed point with probability proportional to
    // the ratio of the actual log density to the log density of the
    // envelope.

    y1 = ds->log_dens(ds, x1, 0, NULL);
    if (log(unif_rand()) <= y1 - eval_segment(nlines, linesegs, x1[which_dim]))
      break;

    // If the proposed point is rejected, but there is no sign of
    // numerical instability, add it to the envelope.  If there is
    // sign of numeric instability, pick a random point in the range
    // of the envelope and add it to the envelope.  This prevents an
    // exceptionally bad envelope from repeatedly proposing and rejecting
    // the same point over and over, which could happen if its points are
    // all deep in sharply decayed tails of a distribution.

    if (R_FINITE(y1) && !degen) {
      add_to_envelope(envelope, x1[which_dim], y1);
    } else if (degen) {
      double U = unif_rand();
      x1[which_dim] =
          U * envelope->absc[0] + (1 - U) * envelope->absc[envelope->nabsc - 1];
      y1 = ds->log_dens(ds, x1, 0, NULL);
      if (R_FINITE(y1)) {
        add_to_envelope(envelope, x1[which_dim], y1);
      }
    }
  }

  // Perform a Metropolis-Hastings transition.  This always accepts
  // if the log density is unimodal.

  double y0 = ds->log_dens(ds, x0, 0, NULL);  // FIXME: could be memoized
  double h0 = eval_segment(nlines, linesegs, x0[which_dim]);
  double h1 = eval_segment(nlines, linesegs, x1[which_dim]);
  Free(linesegs);

  if (log(unif_rand()) > min(0, y1 + min(y0, h0) - y0 - min(y1, h1))) {
    if (rejections != NULL) {
      (*rejections)++;
    }
    return x0[which_dim];
  } else {
    return x1[which_dim];
  }
}

// Create a copy of an envelope.  The caller must call free_envelope()
// on the returned envelope when they are finished with it.

static envelope_t *copy_envelope(envelope_t *env) {
  envelope_t *env2 = Calloc(1, envelope_t);
  env2->nabsc = env->nabsc;
  env2->absc = Calloc(env2->nabsc, double);
  memmove(env2->absc, env->absc, sizeof(double) * env2->nabsc);
  env2->logdens = Calloc(env2->nabsc, double);
  memmove(env2->logdens, env->logdens, sizeof(double) * env2->nabsc);
  return env2;
}

// Free the memory occupied by and envelope by freeing the abscissae
// and log densities, then the envelope itself.

static void free_envelope(envelope_t *envelope) {
  Free(envelope->logdens);
  Free(envelope->absc);
  envelope->logdens = NULL;
  envelope->absc = NULL;
  Free(envelope);
}

// Add a new coordinate pair to an envelope.

static void add_to_envelope(envelope_t *envelope, double new_absc,
                            double new_dens) {
  Rassert(R_FINITE(new_dens));

  // Grow the abscissae and log density array to accommodate the new
  // coordinates.

  envelope->absc = Realloc(envelope->absc, envelope->nabsc + 1, double);
  envelope->logdens = Realloc(envelope->logdens, envelope->nabsc + 1, double);

  // Loop over the abscissae, searching for the first one larger
  // than the new one.

  for (int i = 0; i < envelope->nabsc; i++) {
    if (new_absc < envelope->absc[i]) {
      // We found the next larger abscissa; slide the whole array
      // of abscissae and log densities one place to the right,
      // insert the new one at this index, and return.

      memmove(&envelope->absc[i + 1], &envelope->absc[i],
              sizeof(double) * (envelope->nabsc - i));
      memmove(&envelope->logdens[i + 1], &envelope->logdens[i],
              sizeof(double) * (envelope->nabsc - i));
      envelope->absc[i] = new_absc;
      envelope->logdens[i] = new_dens;
      envelope->nabsc++;
      return;
    }
  }

  // There is no larger abscissa, so the new one must be the new last
  // one, so add it on to the end of the array.

  envelope->absc[envelope->nabsc] = new_absc;
  envelope->logdens[envelope->nabsc] = new_dens;
  envelope->nabsc++;
  return;
}

// Draw a variate from the piecewise exponential distribution
// represented by the nlines-long array of line segments linesegs.  If
// degen is non-NULL, it is set to zero or one indicating whether
// numerical instability was detected (so that the returned point
// should probably not be added to an envelope).

static double sample_piecewise_exponential(int nlines, lineseg_t *linesegs,
                                           int *degen) {
  if (degen != NULL) {
    *degen = 0;
  }

  // Draw a uniform variate.  Iterate through line segments,
  // subtracting the weight of that segment off U until U is less than
  // the current segment's weight.  This chooses a segment randomly
  // with weight proportionally to its weight field.

  double U = unif_rand();
  for (int i = 0; i < nlines; i++) {
    if (U < linesegs[i].weight) {
      // Segment i chosen, draw a variate from that segment.

      double x = sample_chopped_exponential(&linesegs[i]);

      // If the variate is in the first or last thousandth of a
      // percent of the segment, we're probably running into cancellation
      // and should not make the envelope even finer than it already
      // is, so set degen to one.

      double scale = linesegs[i].right - linesegs[i].left;
      if (((x - linesegs[i].left) / scale < 1e-5 ||
           (linesegs[i].right - x) / scale < 1e-5) &&
          degen != NULL) {
        *degen = 1;
      }

      return x;
    } else {
      U -= linesegs[i].weight;
    }
  }
  // Might get here due to roundoff error.  Try again, maybe?
  error("Fatal roundoff error sampling from piecewise exponential.");
  return NAN;
}

// Return a random variate from the exponential density represented
// by the line segment ls.  The log-density of the variate is the line
// with specified slope and intercept conditional on the variate being
// within the bounds of the line segment.  This is used to draw a sample
// from an envelope, once the variate is known to be somewhere in
// an interval.

static double sample_chopped_exponential(lineseg_t *ls) {
  double m = ls->slope;
  double b = ls->icept;
  double l = ls->left;
  double r = ls->right;
  double U = unif_rand();
  Rassert(R_FINITE(m) && R_FINITE(b) && l < r);

  if (!R_FINITE(1 / m)) {
    return (1 - U) * l + U * r;
  }

  double X = (logspace_add(m * l + b + log(1 - U), m * r + b + log(U)) - b) / m;
  Rassert(R_FINITE(X));
  return X;
}

// An envelope_t is a set of (absc,logdens) points, defining an
// envelope by the rules from sec. 2.2 (p. 339) of Gilks & Wild (1992),
// "Adaptive Rejection Sampling for Gibbs Sampling," Applied Statistics
// 41.2:337-348.
//
// The resulting line segments are returned in linesegs_ptr; the
// length of the array is returned in nlines_ptr.  The caller is
// responsible for calling Free() on *linesegs_ptr when they are
// finished with it.

static void envelope_segments(envelope_t *envelope, lineseg_t **linesegs_ptr,
                              int *nlines_ptr) {
  // There must be at least four abscissae.  The segment connecting
  // the first two must slope up, and the segment connecting the last
  // two must slope down.

  Rassert(envelope->nabsc >= 4);
  Rassert(envelope->logdens[0] < envelope->logdens[1]);
  Rassert(envelope->logdens[envelope->nabsc - 1] <
          envelope->logdens[envelope->nabsc - 2]);

  // The first two segments must be special cased.  The first extends
  // the line connecting the first two points from -Inf to the first
  // point.  The second extends the line connecting the second and
  // third to the interval between the first and second, without the
  // possibility of that interval being affected by segments extended
  // from the left of the first.

  // allocate max. possible # of segments
  lineseg_t *linesegs = Calloc(2 * envelope->nabsc, lineseg_t);
  int nlines = 0;

  nlines = 2;

  linesegs[0].left = R_NegInf;
  linesegs[0].right = envelope->absc[0];
  envelope_line(envelope, 0, &linesegs[0].slope, &linesegs[0].icept);

  linesegs[1].left = envelope->absc[0];
  linesegs[1].right = envelope->absc[1];
  envelope_line(envelope, 1, &linesegs[1].slope, &linesegs[1].icept);

  // At iteration i, we process the line segments between abscissae i
  // and i+1 by extending the lines connecting points i-1 and i and
  // points i+1 and i+2.  If one line bounds the other on the interval
  // (absc[i],absc[i+1]), then we just need that segment.  If they
  // intersect over that interval, we need parts of both segments.

  for (int i = 1; i <= envelope->nabsc - 3; i++) {
    double m1, b1, m2, b2;
    double left = envelope->absc[i];
    double right = envelope->absc[i + 1];
    if (left == right) {
      continue;
    }
    envelope_line(envelope, i - 1, &m1, &b1);
    envelope_line(envelope, i + 1, &m2, &b2);
    double midpt = line_intersection(m1, b1, m2, b2);
    if (!R_FINITE(midpt)) {
      midpt = (right + left) / 2;
    }
    Rassert(left < right);
    Rassert((R_FINITE(m1) && R_FINITE(b1)) || (R_FINITE(m2) && R_FINITE(b2)));

    // If the two lines intersect anywhere to the right of the left side of
    // the interval, we need the segment extended from points i-1 and i.

    if (midpt > left || !R_FINITE(m2) || !R_FINITE(b2)) {
      nlines++;
      linesegs[nlines - 1].left = left;
      linesegs[nlines - 1].right = min(midpt, right);
      linesegs[nlines - 1].slope = m1;
      linesegs[nlines - 1].icept = b1;
    }

    // If the two lines intersect anywhere to the left of the right side of
    // the interval, we need the segment extended from points i and i+1.

    if (midpt < right || !R_FINITE(m1) || !R_FINITE(b1)) {
      nlines++;
      linesegs[nlines - 1].left = max(left, midpt);
      linesegs[nlines - 1].right = right;
      linesegs[nlines - 1].slope = m2;
      linesegs[nlines - 1].icept = b2;
    }
  }

  // Just as it was necessary to special case the first two segments,
  // it is necessary to special case the last two, since neither can
  // be affected by lines connecting points to their right.  The second
  // to last line segment is the extension of the line connecting the
  // second and third to last points.  The last line segment is the
  // extension of the line connecting the last two points out to Inf.

  nlines += 2;

  linesegs[nlines - 2].left = envelope->absc[envelope->nabsc - 2];
  linesegs[nlines - 2].right = envelope->absc[envelope->nabsc - 1];
  envelope_line(envelope, envelope->nabsc - 3, &linesegs[nlines - 2].slope,
                &linesegs[nlines - 2].icept);

  linesegs[nlines - 1].left = envelope->absc[envelope->nabsc - 1];
  linesegs[nlines - 1].right = R_PosInf;
  envelope_line(envelope, envelope->nabsc - 2, &linesegs[nlines - 1].slope,
                &linesegs[nlines - 1].icept);

  // Fill in weights.  Each line segment represents the log-density
  // of a piecewise-exponential distribution, with weight proportional
  // to the integral of the exponential of the segment's ordinate over
  // the interval on which the segment is defined.

  // First, compute the normalizing constant, Z,  by adding up the
  // integrals of the exponentiated segment coordinates.  Store the
  // integrals in the weight field of the envelope.

  double Z = R_NegInf;
  for (int j = 0; j < nlines; j++) {
    // Make sure the first and last segments extend to -Inf and +Inf.

    if (j == 0) {
      Rassert(linesegs[j].left == R_NegInf);
      Rassert(linesegs[j].slope > 0);
    } else {
      Rassert(linesegs[j].left == linesegs[j - 1].right);
    }
    if (j == nlines - 1) {
      Rassert(linesegs[j].right == R_PosInf);
      Rassert(linesegs[j].slope < 0);
    }

    // Compute the integral of the line segment analytically.  Avoid
    // explicit exponentials to prevent overflow/underflow, and break
    // apart the cases m>0 and m<0 so that we don't need to take the
    // log of a negative number.

    double m = linesegs[j].slope;
    double b = linesegs[j].icept;
    if (m > 0) {
      linesegs[j].weight = logspace_sub(m * linesegs[j].right + b - log(m),
                                        m * linesegs[j].left + b - log(m));
    } else if (m < 0) {
      linesegs[j].weight = logspace_sub(m * linesegs[j].left + b - log(-m),
                                        m * linesegs[j].right + b - log(-m));
    } else if (m == 0) {
      linesegs[j].weight = b + log(linesegs[j].right - linesegs[j].left);
    } else {
      error("non-finite slope, this should be unreachable");
      linesegs[j].weight = 0;
    }

    Rassert(linesegs[j].weight < R_PosInf);
    Z = logspace_add(Z, linesegs[j].weight);
  }

  // Divide (in log space, subtract) each weight by the normalizing
  // constant so that the weights add up to one.

  for (int j = 0; j < nlines; j++) {
    linesegs[j].weight = exp(linesegs[j].weight - Z);
    Rassert(linesegs[j].weight >= 0);
    Rassert(linesegs[j].weight <= 1.0);
  }

  // Fill in the return values.

  *linesegs_ptr = linesegs;
  *nlines_ptr = nlines;
}

// Returns the line connecting points i and i+1 in envelope.
// If there are numeric issues so that adjacent abscissae have
// the same value, returns a horizontal line.

static void envelope_line(envelope_t *envelope, int i, double *m, double *b) {
  Rassert(i >= 0 && i < envelope->nabsc - 1);
  if (envelope->absc[i + 1] == envelope->absc[i]) {
    *m = 0;
    *b = envelope->logdens[i];
  } else {
    *m = (envelope->logdens[i + 1] - envelope->logdens[i]) /
         (envelope->absc[i + 1] - envelope->absc[i]);
    *b = envelope->logdens[i] - (*m) * envelope->absc[i];
  }
  Rassert(R_FINITE(*m) && R_FINITE(*b));
}

// Returns the value of x at which the two following lines intersect:
//
//   y = m1 * x + b1
//   y = m2 * x + b2

static double line_intersection(double m1, double b1, double m2, double b2) {
  return (b2 - b1) / (m1 - m2);
}

// Given a connected sequence of line segments, find the y value at
// x.  The first line segment should have left==-Inf, the last should
// have right==+Inf, and the others should have their left equal to
// the previous one's right, so that the sequence is connected.

double eval_segment(int nlines, lineseg_t *linesegs, double x) {
  for (int i = 0; i < nlines; i++) {
    if (x <= linesegs[i].right && x >= linesegs[i].left)
      return linesegs[i].slope * x + linesegs[i].icept;
  }
  Rassert(0);
  return NAN;
}

// Expand env so that the first line in the envelope has positive
// slope and the last has negative slope, so that the envelope as a
// whole represents a density that exists.  x0 is a state; which_dim
// specifies which dimension the envelope cuts along.  ds is the density
// being sampled from.

static void expand_envelope(envelope_t *env, dist_t *ds, double *x0,
                            int which_dim) {
  // Copy x0 into x, a scratch state.

  double x[ds->ndim];
  memmove(x, x0, sizeof(double) * ds->ndim);

  // If the first line is positively sloped, use find_decrease to
  // find an ordinate where the density is lower than the second point.
  // If a point outside the support of L is found, use bsearch_function
  // to find an ordinate with finite log-density.  Assumes the log-density
  // goes to -Inf continuously, so this breaks on chi-squared(1) or uniform.

  if (!(env->logdens[1] > env->logdens[0])) {
    x[which_dim] = env->absc[1];
    double w = env->absc[1] - env->absc[0];
    find_decrease(ds, x, which_dim, -w, env->logdens[1], &env->absc[0],
                  &env->logdens[0]);
    if (!R_FINITE(env->logdens[0])) {
      bsearch_function(ds, x, which_dim, env->absc[0], env->absc[1],
                       env->logdens[0], env->logdens[1], &env->absc[0],
                       &env->logdens[0], 1);
    }
  }

  // Repeat the previous procedure, this time ensuring the last line
  // is negatively sloped.

  if (!(env->logdens[env->nabsc - 2] > env->logdens[env->nabsc - 1])) {
    x[which_dim] = env->absc[env->nabsc - 2];
    double w = env->absc[env->nabsc - 1] - env->absc[env->nabsc - 2];
    find_decrease(ds, x, which_dim, w, env->logdens[env->nabsc - 2],
                  &env->absc[env->nabsc - 1], &env->logdens[env->nabsc - 1]);
    if (!R_FINITE(env->logdens[env->nabsc - 1])) {
      bsearch_function(ds, x, which_dim, env->absc[env->nabsc - 2],
                       env->absc[env->nabsc - 1], env->logdens[env->nabsc - 2],
                       env->logdens[env->nabsc - 1], &env->absc[env->nabsc - 1],
                       &env->logdens[env->nabsc - 1], 1);
    }
  }

  Rassert(env->logdens[1] > env->logdens[0]);
  Rassert(env->logdens[env->nabsc - 2] > env->logdens[env->nabsc - 1]);
  Rassert(env->absc[1] > env->absc[0]);
  Rassert(env->absc[env->nabsc - 1] > env->absc[env->nabsc - 2]);
}

// On a cut through x in dimension dim, construct an initial four
// point envelope such that the line through the first two points is
// positively sloped and the line through the second two points is
// negatively sloped.  w is an initial scale to search with.  x itself
// will be one of the two middle points.  The returned envelope should
// be free'd with free_envelope.

static envelope_t *identify_bounds(dist_t *ds, double *x, int dim, double w) {
  double y = ds->log_dens(ds, x, 0, NULL);
  Rassert(R_FINITE(y));

  // Find points to the left and right of x with lower log density.

  double xleft, xright, yleft, yright;
  find_decrease(ds, x, dim, -w, y, &xleft, &yleft);
  find_decrease(ds, x, dim, w, y, &xright, &yright);

  // If the located points do not have finite log density, use
  // bsearch_function to locate points that do.

  if (!R_FINITE(yleft)) {
    bsearch_function(ds, x, dim, xleft, x[dim], yleft, y, &xleft, &yleft, 1);
  }
  if (!R_FINITE(yright)) {
    bsearch_function(ds, x, dim, x[dim], xright, y, yright, &xright, &yright,
                     1);
  }

  // Allocate an envelope and fill in the first and fourth points.

  envelope_t *env = Calloc(1, envelope_t);
  env->nabsc = 4;
  env->absc = Calloc(4, double);
  env->absc[0] = xleft;
  env->absc[3] = xright;
  env->logdens = Calloc(4, double);
  env->logdens[0] = yleft;
  env->logdens[3] = yright;

  // Find a second interior point with bsearch_function.  The
  // decision whether to make x the second or third point based on
  // yright<yleft is somewhat arbitrary; either segment [xlower,x] or
  // [x,xupper] would work.

  double xextra, yextra;
  if (yright < yleft) {
    env->absc[2] = x[dim];
    env->logdens[2] = y;
    bsearch_function(ds, x, dim, xleft, x[dim], yleft, y, &xextra, &yextra, 0);
    env->absc[1] = xextra;
    env->logdens[1] = yextra;
  } else {
    env->absc[1] = x[dim];
    env->logdens[1] = y;
    bsearch_function(ds, x, dim, x[dim], xright, y, yright, &xextra, &yextra,
                     0);
    env->absc[2] = xextra;
    env->logdens[2] = yextra;
  }

  return env;
}

// Searches the space of ds cut through x in dimension dim on the
// interval [lower,upper] for an ordinate with finite log density
// bounded below by the log densities at lower and upper.  If upper_bound
// is set, the log density at the returned ordinate must be bounded
// above by the larger of the log densities at lower and upper as well.
// If known, these can be passed in in the variables ylower and yupper;
// otherwise, those should be set to NAN and they will be computed.
// The found ordinate is returned in q, with its log-density in y.
// This function assumes the log density is continuous.

static void bsearch_function(dist_t *ds, double *x, int dim, double lower,
                             double upper, double ylower, double yupper,
                             double *q, double *y, int upper_bound) {
  // Scratch coordinates

  double x1[ds->ndim];
  memmove(x1, x, sizeof(double) * ds->ndim);

  // Fill in ylower and yupper if necessary.

  if (isnan(ylower)) {
    x1[dim] = lower;
    ylower = ds->log_dens(ds, x1, 0, NULL);
    Rassert(!ISNAN(ylower));
  }

  if (isnan(yupper)) {
    x1[dim] = upper;
    yupper = ds->log_dens(ds, x1, 0, NULL);
    Rassert(!ISNAN(yupper));
  }

  // If the midpoint of lower and upper has a log density larger
  // than the smaller of ylower and yupper, and is either smaller than
  // the larger of the two or upper_bound is unset, then return the
  // midpoint.

  x1[dim] = 0.5 * upper + 0.5 * lower;
  Rassert(x1[dim] != upper && x1[dim] != lower);
  double y1 = ds->log_dens(ds, x1, 0, NULL);
  Rassert(!ISNAN(y1));

  if ((ylower < yupper && ylower < y1 && (!upper_bound || y1 < yupper)) ||
      (ylower > yupper && yupper < y1 && (!upper_bound || y1 < ylower))) {
    *q = x1[dim];
    *y = y1;
    return;
  }

  // Otherwise, either [lower,midpoint] or [midpoint,upper] must
  // contain a suitable point if L is continuous.  Pick the appropriate
  // segment and call bsearch_function recursively.

  if ((ylower > yupper && y1 <= yupper) || (ylower < yupper && y1 >= yupper))
    bsearch_function(ds, x, dim, lower, x1[dim], ylower, yupper, q, y,
                     upper_bound);
  else if ((ylower > yupper && y1 >= ylower) ||
           (ylower < yupper && y1 <= ylower))
    bsearch_function(ds, x, dim, x1[dim], upper, ylower, yupper, q, y,
                     upper_bound);
  else
    error("non-finite bounds in bsearch_function, should be unreachable");
  return;
}

// Searches the log-density surface on a cut through x in dimension dim
// for a point with smaller log-density than x.  y is a pre-computed
// value of the log density at x.  w is an initial increment off x to
// search.  It is repeatedly doubled until a suitable point is found.
// To search to the left, w should be negative; to search to the right,
// it should be positive.  On return, xdec is the ordinate of the found
// point; ydec is the log density at that point, so that if x[dim] were
// set to *xdec, ds->log_dens would compute *ydec.

static void find_decrease(dist_t *ds, double *x, int dim, double w, double y,
                          double *xdec, double *ydec) {
  Rassert(y > R_NegInf);
  Rassert(R_FINITE(w));
  double x1[ds->ndim];
  memmove(x1, x, sizeof(double) * ds->ndim);
  x1[dim] += w;
  double y1 = ds->log_dens(ds, x1, 0, NULL);
  Rassert(!ISNAN(y1));
  if (y1 < y) {
    *xdec = x1[dim];
    *ydec = y1;
    return;
  } else {
    find_decrease(ds, x, dim, 2 * w, y, xdec, ydec);
    return;
  }
}
