% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{make.cone.dist}

\alias{make.cone.dist}

\title{Create a cone distribution object}

\description{Create a cone distribution object as defined by Roberts
  and Rosenthal}

\usage{make.cone.dist(ndim)}

\arguments{
  \item{ndim}{The dimension of the distribution's state space.}
}

\details{
  Defines a distribution object with the following log density:
  \deqn{\pi(x)=e^{-\|x\|}}{pi(x) = exp(-||x||)}
  This is used in Roberts and Rosenthal (2002) to demonstrate
  the deterioration in slice sampler performance as dimensionality
  increases.  It is intended to be passed to
  \code{\link{compare.samplers}}.  Its implementation also serves
  as a simple demonstration of how to define a distribution in C.
}

\value{
  A \code{scdist} object.
}

\seealso{
  \code{\link{make.dist}}
}

\references{
  Roberts, G. O. and Rosenthal, J. S. (2002) \dQuote{The Polar Slice
  Sampler,} Stochastic Models 18(2):257-280.
}
