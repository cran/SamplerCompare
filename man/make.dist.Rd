% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{make.dist}

\alias{make.dist}

\title{Define a probability distribution object}

\description{Defines a probability distribution object for use with
\code{\link{compare.samplers}}.}

\usage{make.dist(ndim, name, name.expression=NULL,
          log.density=NULL, grad.log.density=NULL,
          log.density.and.grad=NULL, initial=NULL,
          mean=NULL, cov=NULL, mean.log.dens=NULL)}

\arguments{
  \item{ndim}{The size of the distribution's state space.}
  \item{name}{A human-readable name for the distribution.}
  \item{name.expression}{A name for the distribution in
    \code{plotmath} notation.  Used in preference to \code{name} in
    plot functions when available.}
  \item{log.density}{A function taking a vector argument that returns
    the log density of the distribution evaluated at that point.}
  \item{grad.log.density}{A function taking a vector argument that
    returns the gradient of the log density of the distribution
    evaluated at that point.}
  \item{log.density.and.grad}{A function taking a vector argument
    and a logical that returns a list with two elements,
    \code{log.density} and \code{grad.log.density}.  The logical
    indicates whether the caller wants the gradient; if not, this
    function may omit the \code{grad.log.density} element in the
    return value.}
  \item{initial}{A function that returns an overdispersed initial
    state for an MCMC simulation of this distribution, used by
    \code{\link{compare.samplers}}.  If unset, uniform draws on a
    unit hypercube are assumed to be acceptable.}
  \item{mean}{A vector specifying the true mean of the distribution.}
  \item{cov}{A matrix specifying the true covariance of the distribution.}
  \item{mean.log.dens}{A scalar specifying the true mean of the log
    density of the distribution.  This will depend on the normalization
    of the log density function.}
}

\details{
  Every distribution must have a name and a dimension.  The log
  density and its gradient are optional; they are used by samplers
  implemented in R.  Samplers implemented in other languages could
  specifically recognize the name of the distribution instead of
  calling back into R, though there is a mechanism for C functions
  to call back.  The mean and covariance do not affect sampling,
  only post-sample diagnostics like autocorrelation time.

  For many distributions, it is easier to compute the log density
  and its gradient at the same time than separately; these will
  generally specify \code{log.density.and.grad} and leave
  \code{log.density} and \code{log.density.and.grad} as \code{NULL}.  The
  returned object will fill those in with calls to
  \code{log.density.and.grad}.  Similarly, if it is simpler to
  compute them separately, \code{log.density.and.grad} will be
  synthesized from \code{log.density} and \code{grad.log.density}
  if necessary.

  \code{mean}, \code{cov}, and \code{mean.log.dens} values are
  intended to be used by diagnostic routines.  \code{mean} and
  \code{mean.log.dens} are currently used by \code{\link{compare.samplers}}
  when estimating autocorrelation times.

  See \code{\link{make.c.dist}} for a way to define distributions
  whose densities are implemented in C instead of R.
}

\value{
  A \code{scdist} object.  It has elements with the same names as the arguments
  to \code{make.dist}.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{make.c.dist}},
  \code{\link{check.dist.gradient}},
  \dQuote{R/C Glue in SamplerCompare} (vignette)
}

\examples{
  # A one dimensional Gamma(3,2) distribution.

  # So that the density does not return NaN outside the support.
  inflog <- function(x) ifelse(x<=0, -Inf, log(x))

  # Define density; unnormalized densities are fine.

  gamma32.log.density <- function(x) (3-1)*inflog(x) - x/2
  gamma32.grad <- function(x) (3-1)/x - 1/2

  # Use make.dist to define the distribution object.

  gamma32.dist <- make.dist(1, 'Gamma32', 'plain("Gamma")(3,2)',
                            log.density=gamma32.log.density,
                            grad.log.density=gamma32.grad,
                            mean=3*2, cov=as.matrix(3*2^2))

  # Make sure the log density and gradient agree at an arbitrary point.

  check.dist.gradient(gamma32.dist, 17)
}
