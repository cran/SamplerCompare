% From SamplerCompare, (c) 2010 Madeleine Thompson

\name{arms.sample}

\alias{arms.sample}

\title{Adaptive Rejection Metropolis Sampler}

\description{Generate a sample from a probability distribution with
  Adaptive Rejection Metropolis Sampling}

\usage{arms.sample(target.dist, x0, sample.size, tuning=1)}

\arguments{
  \item{target.dist}{Target distribution; see \code{\link{make.dist}}.}
  \item{x0}{Numeric vector containing initial state.}
  \item{sample.size}{Sample size requested.}
  \item{tuning}{Scale for initial envelope; see details.}
}

\details{
  \code{arms.sample} implements Adaptive Rejection Metropolis
  Sampling (Gilks, Best, and Tan, 1995).  As described by Gilks et
  al, a user of ARMS must specify an initial envelope roughly
  approximating the target density.  This implementation attempts
  to provide a simpler interface for users by generating an envelope
  automatically.

  To form an initial envelope for coordinate (i), four abscissae
  are needed.  One is \code{x0}.  The sampler tries points
  with abscissae \code{x0[i]-2^k*tuning} and \code{x0[i]+2^k*tuning} for
  whole-numbers \code{k} until points with log densities smaller
  than that at \code{x0} are found, then chooses a fourth point
  from the interior of the two found points.  (Specifically, the
  interval between \code{x0} and the lowest density found point is
  binary-searched until a point with log-density larger than the
  found point is located.)

  This scheme for defining an envelope does not depend on the current
  state in the dimension being sampled.  For discussion of why this
  must be the case, see see Gilks, Neal, Best and Tan (1997).
}

\value{
  A list with elements \code{X}, \code{evals}, and \code{rejections},
  following the calling convention of \code{\link{compare.samplers}}.
  \code{rejections} indicates how many Metropolis-Hastings proposals
  were rejected.
}

\seealso{
  \code{\link{compare.samplers}}
}

\references{
  Gilks, W. R., Best, N. G., and Tan, K. K. C. (1995) \dQuote{Adaptive
  Rejection Metropolis Sampling within Gibbs Sampling,} Applied
  Statistics 44(4):455-472.

  Gilks, W. R., Neal, R. M., Best, N. G., and Tan, K. K. C. (1997)
  \dQuote{Corrigendum: Adaptive Rejection Metropolis Sampling,}
  Applied Statistics 46(2):541-542.
}
