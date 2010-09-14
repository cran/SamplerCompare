% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: cov.match.sample.Rd 1607 2010-09-11 15:45:13Z mthompson $

\name{cov.match.sample}

\alias{cov.match.sample}

\title{Sample with covariance-matching slice sampling}

\description{Generate a sample from a probability distribution with
  the covariance-matching slice sampling method.}

\usage{cov.match.sample(target.dist, x0, sample.size, tuning=1,
                 theta=1, limit=length(x0)*100)}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Initial coordinates.}
  \item{sample.size}{Sample size to draw.}
  \item{tuning}{A tuning parameter; corresponds to \eqn{\sigma_c}
    in sec. 4 of Thompson and Neal (2010).}
  \item{theta}{A factor to scale the crumb standard deviation in
    every direction after a proposal is rejected.  So, after \eqn{k}
    proposals, crumbs have standard deviation \code{theta^k*tuning}
    in directions orthogonal to all the proposal gradients.}
  \item{limit}{A limit on the number of log-density evaluations per
    observation before sampling is aborted.}
}

\details{
  This function implements the covariance-matching method of slice
  sampling, as described by Thompson and Neal (2010).  It can be passed
  to \code{\link{compare.samplers}} in the \code{samplers} list
  argument.
}

\value{
  A list with elements \code{X}, \code{evals}, \code{grads}, and
  \code{adapt.rate}.  \code{adapt.rate} indicates the fraction of
  crumb draws that resulted in adaptation.  This sampler follows
  the calling convention of \code{\link{compare.samplers}}.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{shrinking.rank.sample}}
}

\references{
  Thompson, M. B. and Neal, R. M. (2010). Covariance-adaptive slice
  sampling. Technical Report TR-1002, Dept. of Statistics, University
  of Toronto.
}
