% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{hyperrectangle.sample}

\alias{hyperrectangle.sample}
\alias{nograd.hyperrectangle.sample}

\title{Multivariate slice samplers}

\description{Generate a sample from a probability distribution with
  a slice sampler taking multivariate steps.}

\usage{hyperrectangle.sample(target.dist, x0, sample.size, tuning=1,
                 use.gradient=TRUE, limit=length(x0)*100)
nograd.hyperrectangle.sample(...)
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Sample size requested.}
  \item{tuning}{Initial edge length of hyperrectangle.}
  \item{use.gradient}{A logical indicating whether the sampler
    should use the gradient when shrinking the box.}
  \item{limit}{A limit on the number of log-density evaluations per
    observation before sampling is aborted.}
  \item{...}{\code{nograd.hyperrectangle.sample} takes the same arguments
    as \code{hyperrectangle.sample}, except \code{use.gradient}.}
}

\details{
  \code{hyperrectangle.sample} implements multivariate slice sampling
  with hyperrectangles as described in Neal (sec. 5.1, 2003).
  
  If \code{use.gradient} is set, when a proposal is rejected, the
  gradient at the rejected proposal is used to choose a direction
  to shrink the box.  Neal suggested shrinking in the direction the
  gradient was largest, but this implementation shrinks in the
  direction that the gradient times the box length is largest to
  better handle poorly scaled distributions.

  If \code{use.gradient} is not set, the gradient is not computed
  and the box is shrunk in all directions after every rejected
  proposal.  Calling \code{nograd.hyperrectangle.sample} is equivalent
  to calling \code{hyperrectangle.sample} with \code{use.gradient=FALSE};
  the extra name is provided for convenience when using either of
  these functions with \code{\link{compare.samplers}}.
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{grads}.  
  This sampler follows the calling convention of
  \code{\link{compare.samplers}}.
}

\references{
  Neal, Radford M. (2003), \dQuote{Slice Sampling,} The Annals of
  Statistics 31(3):705-767.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{nonadaptive.crumb.sample}},
  \code{\link{interval.slice.sample}}
}
