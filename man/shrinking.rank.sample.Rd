% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{shrinking.rank.sample}

\alias{shrinking.rank.sample}

\title{Sample with shrinking-rank slice sampling}

\description{Generate a sample from a probability distribution with
  the shrinking-rank slice sampling method.}

\usage{shrinking.rank.sample(target.dist, x0, sample.size, tuning=1,
         downscale=0.95, min.dimension=1)}

\details{
  \code{shrinking.rank.slice.sample}  implements the shrinking-rank
  method of slice sampling, as described by Thompson and Neal (2010).  It
  can be passed to \code{\link{compare.samplers}} in the \code{samplers}
  list argument.
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Requested sample size.}
  \item{tuning}{A tuning parameter; corresponds to \eqn{\sigma_c}
    in sec. 5 of Thompson and Neal (2010).}
  \item{downscale}{Factor to reduce crumb standard deviation by when a
    proposal is rejected.}
  \item{min.dimension}{The minimum dimension to sample crumbs from.}
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{grads},
  following the calling convention of \code{\link{compare.samplers}}.
}

\seealso{
  \code{\link{compare.samplers}}
  \code{\link{cov.match.sample}}
}

\references{
  Thompson, M. B. and Neal, R. M. (2010). Covariance-adaptive slice
  sampling. Technical Report TR-1002, Dept. of Statistics, University
  of Toronto.
}
