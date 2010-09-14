% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: check.dist.gradient.Rd 1484 2010-08-24 15:39:22Z mthompson $

\name{check.dist.gradient}

\alias{check.dist.gradient}

\title{Test a gradient function}

\description{Evaluates the gradient function of a distribution
  object and compares it to a numeric gradient computed from the log
  density function.}

\usage{check.dist.gradient(ds, x, h=1e-7)}

\arguments{
  \item{ds}{A distribution object with defined log density and
    gradient functions.}
  \item{x}{A point inside the support of \code{ds}.}
  \item{h}{An offset from \code{x} at which to evaluate the log
    density when computing numeric derivatives.}
}

\details{
  \code{check.dist.gradient} computes the numeric derivative of
  \code{ds$log.density} at \code{x} in each of its coordinates and
  compares this to the value returned by \code{ds$grad.log.density}.
  If the relative error is greater than \eqn{0.001}, an error is
  reported.

  This function can be used when defining a distribution to ensure
  that the gradient function is implemented correctly.
}

\examples{
  check.dist.gradient(N2weakcor.dist, runif(2))
}

\seealso{
  \code{\link{make.dist}}
}
