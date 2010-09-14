% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: adaptive.metropolis.sample.Rd 1490 2010-08-26 02:50:46Z mthompson $

\name{adaptive.metropolis.sample}

\alias{adaptive.metropolis.sample}

\title{Adaptive Metropolis}

\description{Generate a sample from a probability distribution with
  the Adaptive Metropolis algorithm}

\usage{adaptive.metropolis.sample(target.dist, x0, sample.size,
                           tuning=0.1, beta=0.05, burn.in=0.2)
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Requested sample size.}
  \item{tuning}{Standard deviation of first component of proposal
    distribution}
  \item{beta}{Weight of first component of proposal distribution}
  \item{burn.in}{Stop adaptation after this fraction of the chain.
    Set this to 1.0 to obtain the behavior described by Roberts and
    Rosenthal (2009).}
}

\details{
  This function implements the Adaptive Metropolis algorithm as
  described by Roberts and Rosenthal (2009).  Proposals are a mixture
  of a spherical Gaussian with standard deviation equal to
  \code{tuning/sqrt(target.dist$ndim)} (with weight \code{beta})
  and a Gaussian with covariance equal to the sample covariance of
  the already-computed observations scaled by
  \code{2.38^2/target.dist$ndim} (with weight \code{1-beta}).
  The resulting Markov chain is not strictly stationary with the
  target distribution for the burn-in period of the chain, but
  is ergodic.
}

\value{
  A list containing the elements \code{X}, \code{evals},
  \code{reject.rate}, and \code{sample.cov}.  This sampler follows
  the calling convention of \code{\link{compare.samplers}}.
  \code{reject.rate} contains the fraction of proposals that were
  rejected.  \code{sample.cov} is the most recent sample covariance
  used to update the proposal distribution.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{multivariate.metropolis.sample}}
}

\references{
  Roberts, G. O. and Rosenthal, J. S. (2009), \dQuote{Examples of
  Adaptive MCMC,} Journal of Computational and Graphical Statistics
  18(2):349-367.
}
