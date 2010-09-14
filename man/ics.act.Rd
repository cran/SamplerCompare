% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: ics.act.Rd 1607 2010-09-11 15:45:13Z mthompson $

\name{ics.act}

\alias{ics.act}

\title{Compute the autocorrelation time of a chain}

\description{Computes the autocorrelation time of an MCMC chain
  using the initial convex sequence.}

\usage{ics.act(X, true.mean=colMeans(as.matrix(X)))}

\arguments{
  \item{X}{A matrix or vector containing the states of a stationary
    Markov chain.  Each row is a single state.}
  \item{true.mean}{A vector containing the true mean of \code{X}.
    It should have as many elements as \code{X} has columns.  The
    default is to use the sample mean; this can lead the procedure
    to underestimate the ACT if the simulation has not converged.}
}

\details{
  \code{ics.act} computes the autocorrelation time of chain represented
  by \code{X} using the initial convex sequence (Geyer 1992).  It
  is roughly equivalent to the maximum of \code{initseq(xx)$var.con
  / initseq(xx)$gamma0} where \code{xx} ranges over the columns of
  \code{X}, but properly handles the case where the true mean of
  the target distribution is known.

  Callers may want to remove a burn-in period from a sample before
  passing it to \code{ics.act}.

  This is not intended for regular use; it is included as a way to
  cross-check the results of \code{ar.act}.  It depends on the
  \pkg{mcmc} package.
}

\value{
  The estimated autocorrelation time of \code{X}.
}

\note{
  This function bypasses the interface of \pkg{mcmc}, and may break
  if the implementation of \code{initseq} changes.
}

\seealso{
  \code{\link{ar.act}},
  \code{mcmc::initseq}
}

\references{
  Geyer, C. J. (1992) \dQuote{Practical Markov Chain Monte Carlo,}
  Statistical Science 7(4):473-511.
}
