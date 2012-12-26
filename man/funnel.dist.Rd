% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{funnel.dist}

\alias{funnel.dist}

\title{Funnel distribution object}

\description{A distribution object for Radford Neal's funnel distribution}

\details{
  \code{funnel.dist} represents the funnel distribution described by Neal
  (2003, p. 732).  It is a ten-dimensional distribution on the reals, with:
  \deqn{v \sim N(0,3^2),\; x[k] \sim N(0,e^v) \;\mathrm{for}\; k=1,\ldots,9}{%
        v~N(0,3^2), x[k]~N(0,exp(v)) for k=1,...,9}
  The state space is \eqn{(v,x[1],x[2],...,x[9])}.  The name comes
  from the funnel-shaped two dimensional marginal distributions
  \eqn{(v,x[k])}.

  This object is intended as a demonstration to be passed to
  \code{\link{compare.samplers}}.
}

\references{
  Neal, Radford M. (2003), \dQuote{Slice Sampling,} The Annals of
  Statistics 31(3):705-767.
}

\seealso{
  \code{\link{make.dist}},
  \code{\link{compare.samplers}}
}
