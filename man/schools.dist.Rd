% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{schools.dist}

\alias{schools.dist}

\title{Eight schools distribution object}

\description{A distribution object for the eight-schools distribution}

\details{
  This object represents the distribution of \dQuote{eight schools,}
  a ten-dimensional multilevel model from Gelman et al (2004).  The
  first and second parameters are mean and log-variance
  hyperparameters, and the third through tenth are group-level
  means.

  This object is intended as a demonstration to be passed to
  \code{\link{compare.samplers}}.
}

\seealso{
  \code{\link{make.dist}},
  \code{\link{compare.samplers}}
}

\references{
  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (2004).
  Bayesian Data Analysis, Second Edition. Chapman and Hall/CRC.
  pp. 138-145.
}
