% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: chud.Rd 1490 2010-08-26 02:50:46Z mthompson $

\name{chud}

\alias{chud}
\alias{chdd}

\title{Cholesky Update/Downdate}

\description{Rank-one updates of Cholesky factors}

\usage{chud(R,x)
chdd(R,x)}

\arguments{
  \item{R}{an upper-triangular matrix}
  \item{x}{a vector}
}

\value{
  An updated version of \code{R}.
}

\details{
  \code{chud} computes Q such that: \deqn{Q^T Q = R^T R + x x^T}
  \code{chdd} computes Q such that: \deqn{Q^T Q = R^T R - x x^T}
  \code{chdd} reports an error if \eqn{R^T R - x x^T} is not positive
  definite.  The two functions use LINPACK's \code{dchud} and
  \code{dchdd} routines respectively, two of the few routines from
  LINPACK without analogues in LAPACK.
}

\seealso{
  \code{\link{chol}}
}

\references{
  Dongarra, J. J., Moler, C. B., Bunch, J. R., Stewart, G. W. (1979)
  LINPACK User's Guide.
}
