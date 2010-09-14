% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: multivariate.metropolis.sample.Rd 1490 2010-08-26 02:50:46Z mthompson $

\name{multivariate.metropolis.sample}

\alias{multivariate.metropolis.sample}
\alias{univar.metropolis.sample}

\title{Metropolis samplers}

\description{Generate a sample from a probability distribution with
  the Metropolis algorithm.}

\usage{multivariate.metropolis.sample(target.dist, x0, sample.size, tuning=1)
univar.metropolis.sample(target.dist, x0, sample.size, tuning=1)
}

\details{
  These two functions implement variants of the Metropolis algorithm
  for sampling a target distribution, following the interface
  used by \code{\link{compare.samplers}}.
  \code{multivariate.metropolis.sample} uses spherically symmetric
  Gaussian proposals with marginal standard deviation equal to the
  \code{tuning} parameter.  \code{univar.metropolis.sample} updates
  each coordinate in sequence using univariate Gaussian proposals
  with standard deviation equal to the \code{tuning} parameter.
  
  So that these two functions are roughly comparable, with a
  \eqn{p}-dimensional target distribution,
  \code{multivariate.metropolis.sample} performs \eqn{p} accept-reject
  steps each time between observations, so that both functions
  evaluate the log density a number of times roughly equal to \eqn{p}
  times the sample size.  While there are often efficiency optimizations
  possible when only one coordinate is updated,
  \code{univar.metropolis.sample} does not support these; these two
  samplers are included for comparison rather than for practical use.
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Sample size requested.}
  \item{tuning}{Proposal standard deviation}
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{reject.rate}.
  See \code{\link{compare.samplers}} for more information on \code{X} and
  \code{evals}.  \code{reject.rate} is the fraction of proposals
  not accepted.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{adaptive.metropolis.sample}}
}
