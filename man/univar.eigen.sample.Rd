% From SamplerCompare, (c) 2011 Madeleine Thompson

\name{univar.eigen.sample}

\alias{univar.eigen.sample}
\alias{cheat.univar.eigen.sample}

\title{Eigendecomposition-based slice samplers}

\description{Generate a sample from a probability distribution with
  slice sampling with univariate steps along eigenvectors.}

\usage{univar.eigen.sample(target.dist, x0, sample.size, tuning=1,
                    steps.out=100, cheat=FALSE)
cheat.univar.eigen.sample(target.dist, x0, sample.size, tuning=1,
                    steps.out=100)
}

\details{
  These two functions implement slice sampling with univariate steps
  along estimated eigenvectors.  Thompson (2011, ch. 3) has details
  on the algorithms.  The functions follow the interface used by
  \code{\link{compare.samplers}}.  Calling \code{cheat.univar.eigen.sample}
  is equivalent to calling \code{univar.eigen.sample} with
  \code{cheat=TRUE}; it is provided as a convenience so that it can
  be passed directly to \code{compare.samplers}.
}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Sample size requested.}
  \item{tuning}{Initial slice approximation length.}
  \item{steps.out}{Maximum number of iterations the stepping out
    algorithm should run when choosing an initial slice approximation.
    Set to NULL to refrain from stepping out.}
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
  \code{\link{oblique.hyperrect.sample}}
}

\references{
  Thompson, M. B. (2011), Slice Sampling with Multivariate Steps
  (forthcoming).
}
