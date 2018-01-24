% From SamplerCompare, (c) 2010-2011 Madeleine Thompson

\name{compounded.sampler}

\alias{compounded.sampler}

\title{Build a sampler from transition functions}

\description{Defines a probability distribution object for use with
\code{\link{compare.samplers}}.}

\usage{compounded.sampler(step.functions, name, name.expr=NULL)}

\arguments{
  \item{step.functions}{A list of transition functions; see details.}
  \item{name}{A character string naming the sampler.}
  \item{name.expr}{A character string naming the sampler in
    \code{plotmath} notation.}
}

\details{
  \code{compounded.sampler} builds an MCMC sampler following the
  conventions of \code{\link{compare.samplers}} from a list of
  transition functions.  The returned sampler has four arguments:
  \var{target.dist}, \var{x0}, \var{sample.size}, and \var{limit}.
  Further arguments, including the standard argument \var{tuning},
  are passed to every transition function.  The first three arguments
  transition functions are passed are \var{target.dist}, a vector
  state \var{x} to transition from, and the log density at that
  state, \var{y}.  They should return a list containing four elements:
  \var{x}, \var{y}, \var{evals}, and \var{grads}.  \var{x} is the
  state transitioned to, \var{y} is the log density at that state,
  and \var{evals} and \var{grads} are the number of log density and
  gradient calls made in that transition.

  Each MCMC iteration, the first transition function is called with
  the current state of the chain.  The state it returns is passed
  to the second transition function, whose returned state is passed
  to the third, and so on.  The state returned by the final transition
  function is taken to be the state of the chain as a whole at the
  end of the iteration.

  This way, transition functions that provide complementary features,
  such as fast mixing in different coordinates, can be combined
  without modifying their internal structure.  The \code{transition_fn}
  interface provides a similar mechanism for samplers implemented
  in C.  It is documented in the vignette, \dQuote{R/C Glue in
  SamplerCompare}.
}

\value{
  A sampler function.
}

\seealso{
  \code{\link{compare.samplers}}
}
