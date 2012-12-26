% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{make.gaussian}

\alias{make.gaussian}
\alias{N2weakcor.dist}
\alias{N4poscor.dist}
\alias{N4negcor.dist}

\title{Gaussian distribution objects}

\description{Gaussian distribution objects}

\usage{make.gaussian(mean, sigma=NULL, rho=NULL)
N2weakcor.dist
N4poscor.dist
N4negcor.dist
}

\arguments{
  \item{mean}{The mean of the distribution as a numeric vector;
    implicitly specifies the dimension.}
  \item{sigma}{The covariance of the distribution.}
  \item{rho}{The marginal correlations between parameters.}
}

\details{
  \code{make.gaussian} returns a distribution object representing
  a multivariate normal distribution.  If \code{sigma} is specified,
  that is taken to be its covariance.  Otherwise, if \code{rho} is
  specified, the covariance is taken to be a matrix with ones on
  the diagonal and \code{rho} on the off-diagonal elements.  To
  preserve positive definiteness, \code{rho} must be between
  \code{-1/(length(mean)-1)} and 1.

  \code{N2weakcor.dist}, \code{N4poscor.dist}, and \code{N4negcor.dist}
  are predefined distributions generated with \code{make.gaussian}.
  They are intended to be used as test cases with
  \code{\link{compare.samplers}}.  The examples below show how they
  are defined.  \code{N2weakcor.dist} is a weakly positively
  correlated two-dimensional Gaussian.  \code{N4poscor.dist} is a
  highly positively correlated four-dimensional Gaussian.
  \code{N4negcor.dist} is a highly negatively correlated four-dimensional
  Gaussian.  \code{N4poscor.dist} and \code{N4negcor.dist} are
  similarly conditioned, but \code{N4poscor.dist} has one large
  eigenvalue and three small ones, while \code{N4negcor.dist} has
  one small eigenvalue and three large ones.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link{make.dist}}
}

\examples{
  N2weakcor.dist <- make.gaussian(c(0,0), rho=0.8)
  N4poscor.dist <- make.gaussian(c(1,2,3,4), rho=0.999)
  N4negcor.dist <- make.gaussian(c(1,2,3,4), rho=-0.3329)
}
