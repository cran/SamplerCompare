% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: make.mv.gamma.dist.Rd 1490 2010-08-26 02:50:46Z mthompson $

\name{make.mv.gamma.dist}

\alias{make.mv.gamma.dist}

\title{Create a distribution object for a set of uncorrelated
  Gamma distributions}

\description{Create a distribution object for a set of uncorrelated
  Gamma distributions}

\usage{make.mv.gamma.dist(shape, scale=rep(1, length(shape)))}

\arguments{
  \item{shape}{A vector of shape parameters.}
  \item{scale}{A vector of scale parameters.  Must have the same
    length as \code{shape}}
}

\details{
  Defines a distribution object for a multivariate distribution
  where each marginal density is Gamma and uncorrelated with the
  other coordinates.  The log density is therefore equivalent to
  \code{sum(dgamma(x, shape, scale=scale, log=TRUE))}.

  This is included as a test case for comparing how MCMC methods
  perform on asymmetric distributions.
}

\value{
  A \code{dist} object.  For convenience, the \code{shape} and \code{scale}
  elements are filled in with the parameters passed to
  \code{make.mv.gamma.dist}.
}

\seealso{
  \code{\link{make.gaussian}},
  \code{\link{compare.samplers}}
}
