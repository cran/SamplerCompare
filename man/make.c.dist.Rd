% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: make.c.dist.Rd 3139 2012-01-28 19:24:19Z mthompson $

\name{make.c.dist}

\alias{make.c.dist}

\title{Define a probability distribution object with C log-density}

\description{Defines a probability distribution object for use with
\code{\link{compare.samplers}} with log-density implemented in C.}

\usage{make.c.dist(ndim, name, c.log.density, c.context = NULL,
             name.expression = NULL, mean = NULL, cov = NULL)}

\arguments{
  \item{ndim}{The size of the distribution's state space.}
  \item{name}{A human-readable name for the distribution.}
  \item{c.log.density}{A C function returning the log-density and
    gradient of the target distribution.}
  \item{c.context}{An opaque object passed to \code{c.log.density}}
  \item{name.expression}{A name for the distribution in
    \code{plotmath} notation.  Used in preference to \code{name} in
    plot functions when available.}
  \item{mean}{A vector specifying the true mean of the distribution.}
  \item{cov}{A matrix specifying the true covariance of the distribution.}
}

\details{
  See \code{\link{make.dist}} for discussion of \code{ndim},
  \code{name}, \code{name.expression}, \code{mean}, and \code{cov}.

  \code{c.log.density} is a string containing the symbol name of a
  C function that computes the log density and log density gradient
  of the target distribution.  It has the type \code{log_density_t},
  defined in \code{SamplerCompare.h} as:\preformatted{typedef double log_density_t(dist_t *ds, double *x,
                             int compute_grad, double *grad);
typedef struct {
  log_density_t *log_dens;
  SEXP context;
  int ndim;
} dist_t;}
  The \code{ds} structure defines the distribution, where the
  \code{log_dens} element is a pointer to the function named by
  \code{c.log.density}, the \code{context} element is a \code{SEXP}
  containing the \code{c.context} parameter, and \code{ndim} is the
  \code{ndim} parameter to \code{make.c.dist}.
  
  The \code{x} parameter is an \code{ndim}-long array of doubles
  containing the location at which to evaluate the log-density,
  which the \code{log_density_t} should return.

  If \code{compute_grad} is nonzero, the function should compute
  the gradient of the log density and store it in the double array
  pointed to by \code{grad}.  If for some reason it cannot do this,
  it should call the R-internal \code{error} function to report an
  error to the user.  If the implementor does not plan to sample
  from the distribution with a method that computes gradients, this
  can reduce implementation effort.

  The details of this interface are described in greater detail in
  \dQuote{R/C Glue in SamplerCompare}.
}

\value{
  A \code{dist} object.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{make.dist}},
  \dQuote{R/C Glue in SamplerCompare} (vignette)
}
