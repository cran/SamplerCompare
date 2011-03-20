% From SamplerCompare, (c) 2011 Madeleine Thompson

\name{oblique.hyperrect.sample}

\alias{oblique.hyperrect.sample}
\alias{cheat.oblique.hyperrect.sample}

\title{Eigendecomposition-based hyperrectangle method}

\description{Generate a sample from a probability distribution with
  the hyperrectangle method with slice approximation axes oriented
  along eigenvectors.}

\usage{oblique.hyperrect.sample(target.dist, x0, sample.size, tuning=1,
                    edge.scale=5, cheat=FALSE)
cheat.oblique.hyperrect.sample(target.dist, x0, sample.size, tuning=1)
}

\details{
  These two functions implement the hyperrectangle method (Neal,
  2003, sec. 5.1) with the hyperrectangle oriented along estimates
  of the eigenvectors of the target distribution's covariance, as
  described by Thompson (2011, ch. 3).  The functions follow the
  interface used by \code{\link{compare.samplers}}.  Calling
  \code{cheat.oblique.hyperrect.sample} is equivalent to calling
  \code{oblique.hyperrect.sample} with \code{cheat=TRUE}; it is
  provided as a convenience so that it can be passed directly to
  \code{compare.samplers}.
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Sample size requested.}
  \item{tuning}{Scale of initial/fallback hypercube edge; \var{w}
    in Thompson (2011, ch. 3).}
  \item{edge.scale}{The initial slice approximation has edges of length
    equal to the square root of the corresponding eigenvalue times
    this factor.}
  \item{cheat}{Set to true to use the covariance from \var{target.dist}
    instead of estimating it.  This is not possible on real problems
    but can be useful for debugging.}
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{grads}.
  See \code{\link{compare.samplers}} for more information.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{univar.eigen.sample}}
}

\references{
  Neal, Radford M. (2003), \dQuote{Slice Sampling,} The Annals of
  Statistics 31(3):705-767.

  Thompson, M. B. (2011), Slice Sampling with Multivariate Steps
  (forthcoming).
}
