% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: make.multimodal.dist.Rd 1486 2010-08-25 21:47:38Z mthompson $

\name{make.multimodal.dist}

\alias{make.multimodal.dist}

\title{Create a distribution object for a random mixture of Gaussians}

\description{Create a distribution object for a random mixture of Gaussians}

\usage{make.multimodal.dist(nmodes, ndim, cube.size)}

\arguments{
  \item{nmodes}{The number of components in the mixture model.}
  \item{ndim}{The dimension of the model.}
  \item{cube.size}{The edge length of the hypercube in which the
    modes are distributed.}
}

\details{
  Defines a distribution object for a mixture of random Gaussians.
  The means of the the \code{nmodes} Gaussians are randomly distributed
  over an \code{ndim}-dimensional hypercube with one corner at the
  origin and the opposite \code{cube.size} away in each positive
  direction.  The same random seed is temporarily set when drawing
  modes, so every time this function is called with the same parameters,
  the resulting distribution is the same.
  
  This is included as a test case for comparing how MCMC methods
  perform on multimodal distributions.
}

\value{
  A \code{dist} object.  For convenience, the \code{modes} element
  is a matrix containing the modes as rows.
}

\seealso{
  \code{\link{make.dist}},
  \code{\link{compare.samplers}}
}
