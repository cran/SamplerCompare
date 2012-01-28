% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: ar.act.Rd 3139 2012-01-28 19:24:19Z mthompson $

\name{ar.act}

\alias{ar.act}

\title{Compute the autocorrelation time of a chain}

\description{Computes the autocorrelation time of an MCMC chain
  using an AR model with order chosen by AIC.}

\usage{ar.act(Y, true.mean=NULL)}

\arguments{
  \item{Y}{A matrix or vector containing the states of a stationary
    Markov chain.  If a matrix, each row is a single state.}
  \item{true.mean}{A vector containing the true mean of \code{Y}.
    It should be either NULL or have as many elements as \code{Y} has
    columns.  If NULL, the sample mean of \code{Y} is used.}
}

\details{
  This function fits an AR(\eqn{p}) model to each component of the
  chain with states \var{Y} using the Yule-Walker method to
  estimate the coefficients and AIC to estimate \eqn{p}.  Let
  \eqn{\pi}{pi} be the vector of estimated AR coefficients
  for column \eqn{i}, and let \eqn{\rho}{rho} be the sample autocorrelation
  function to lag \eqn{p}.  Then, the autocorrelation
  time of the component is estimated as:
    \deqn{\tau_i = \frac{1-\pi^T \rho}{(1-\sum\pi)^2}
      }{tau = (1-pi^T rho)/(1-sum(pi))^2}
  For more discussion of this formula and its associated confidence
  intervals, see Thompson (2010).

  The returned autocorrelation time (and associated confidence
  interval) are the maxima over the columns of \var{Y}.

  Callers may want to remove a burn-in period from a sample before
  passing it to \code{ar.act}.
}

\value{
  A list with elements:
  \itemize{
    \item \var{act}: the estimated autocorrelation time of the
      slowest-mixing column of \var{Y}.
    \item \var{se}: the standard error of \var{act}.
    \item \var{act.025}, \var{act.975}: a nominal 95\% confidence
      interval for \var{act}.  Since the interval is asymmetric about
      \var{act}, the standard error is not sufficient to generate
      these.
    \item \var{order}: The order of the AR model selected (\eqn{p}).
  }
}

\seealso{
  \code{\link{compare.samplers}},
  \code{ar.yw},
  \code{CODA::effectiveSize},
  \code{CODA::spectrum0.ar}
}

\references{
  Thompson, M. B. (2010).  Graphical comparison of MCMC performance
  (forthcoming).
}
