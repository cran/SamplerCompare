% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: nonadaptive.crumb.sample.Rd 1490 2010-08-26 02:50:46Z mthompson $

\name{nonadaptive.crumb.sample}

\alias{nonadaptive.crumb.sample}

\title{Sample with nonadaptive-crumb slice sampling}

\description{Generate a sample from a probability distribution with
  the nonadaptive-crumb slice sampling method.}

\usage{nonadaptive.crumb.sample(target.dist, x0, sample.size,
                        tuning=1, downscale=0.95)}

\details{
  This function implements slice sampling with nonadaptive crumbs.
  Crumbs are Gaussian with spherical covariance starting at
  \code{tuning}, decreasing by \code{downscale} each time a proposal
  is rejected.  More information can be found in sec. 5.2 of Neal
  (2003).  This function can be passed to \code{\link{compare.samplers}}
  in the \code{samplers} list argument.
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Requested sample size.}
  \item{tuning}{Initial crumb standard deviation.}
  \item{downscale}{Factor to reduce crumb standard deviation by when a
    proposal is rejected.}
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{grads},
  following the calling convention of \code{\link{compare.samplers}}.
}

\seealso{
  \code{\link{shrinking.rank.sample}},
  \code{\link{compare.samplers}}
}

\references{
  Neal, Radford M. (2003), \dQuote{Slice Sampling,} The Annals of
  Statistics 31(3):705-767.
}
