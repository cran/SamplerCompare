% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: comparison.plot.Rd 1607 2010-09-11 15:45:13Z mthompson $

\name{comparison.plot}

\alias{comparison.plot}

\title{Plot the results of compare.samplers}

\description{Generates a plot of representing results from
  \code{\link{compare.samplers}}.}

\usage{comparison.plot(RS, xlab=NULL, ylab=NULL, ...)}

\arguments{
  \item{RS}{A data frame in the form returned by
    \code{\link{compare.samplers}}}
  \item{xlab, ylab, ...}{Options to be passed to
    \code{\link[ggplot2]{qplot}}.}
}

\details{
  This function generates a grid of subplots, where each column of
  plots represents a sampler and each row represents a distribution.
  The horizontal axis in each subplot represents the tuning parameter
  passed as \code{tuning} to \code{\link{compare.samplers}}, and
  the vertical axis represents the product \code{RS$evals * RS$act},
  the number of log density evaluations per independent sample
  required for that distribution when simulated by that sampler
  with that tuning parameter.  95\% confidence intervals, covering
  the range \code{[RS$evals * RS$act.025, RS$evals * RS$act.975]},
  are represented by vertical bars.  If \code{RS$evals} or \code{RS$act}
  is missing or infinite, a question mark is plotted instead of the
  default plot character.

  \code{comparison.plot} returns a \pkg{ggplot2} plot object.  If it
  is called non-interactively, one must call \code{\link{print}}
  on the returned object for a plot to be displayed.  To superimpose
  other figures of merit on the plot, one can add \code{geom_*}
  objects to the returned plot object before calling \code{\link{print}}.

  For more discussion of this type of plot, see Thompson (2010).
}

\note{
  This is the only function in SamplerCompare that uses the \pkg{ggplot2}
  package, so it is loaded explicitly by \code{comparison.plot}
  instead of being listed as a package dependency.  This way, compute
  servers calling \code{\link{compare.samplers}} do not need to
  have \pkg{ggplot2} installed on them.
}

\value{
  A \pkg{ggplot2} plot object.
}

\seealso{
  \code{\link{compare.samplers}},
  \code{\link[ggplot2]{qplot}},
  \dQuote{Introduction to SamplerCompare} (vignette)
}

\references{
  Thompson, M. B. (2010). Graphical comparison of MCMC performance.
  In preparation.
}