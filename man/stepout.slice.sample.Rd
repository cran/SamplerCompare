% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: stepout.slice.sample.Rd 1490 2010-08-26 02:50:46Z mthompson $

\name{stepout.slice.sample}

\alias{stepout.slice.sample}
\alias{interval.slice.sample}

\title{Univariate slice samplers}

\description{Generate a sample from a probability distribution with
  a slice sampler.}

\usage{stepout.slice.sample(target.dist, x0, sample.size, tuning=1,
                     step.out=TRUE, limit=length(x0)*100)
interval.slice.sample(...)
}

\details{
  \code{stepout.slice.sample} implements univariate slice sampling
  with stepping out as described in sec. 4 of Neal (2003).  If
  \code{step.out=FALSE} or \code{interval.slice.sample} is called
  instead, no stepping out is performed; the wrapper function
  \code{interval.slice.sample} is provided for convenience when
  calling \code{\link{compare.samplers}}.

  If \code{target.dist} is a multivariate distribution, each step
  of the Markov chain updates each coordinate once in sequence.
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Requested sample size.}
  \item{tuning}{Initial interval length for slice.}
  \item{step.out}{Flag indicating whether to expand the initial
    interval before proposing a new coordinate.}
  \item{limit}{A limit on the number of log-density evaluations per
    observation before sampling is aborted.}
  \item{...}{\code{interval.slice.sample} takes the same arguments
    as \code{stepout.slice.sample}, except \code{step.out}.}
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{grads},
  following the calling convention of \code{\link{compare.samplers}}.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{hyperrectangle.sample}}
}

\references{
  Neal, Radford M. (2003), \dQuote{Slice Sampling,} The Annals of
  Statistics 31(3):705-767.
}
