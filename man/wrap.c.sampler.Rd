% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: wrap.c.sampler.Rd 3139 2012-01-28 19:24:19Z mthompson $

\name{wrap.c.sampler}

\alias{wrap.c.sampler}

\title{Create an R stub function for a sampler implemented in C}

\description{Create an R stub function for a sampler implemented in C.}

\usage{wrap.c.sampler(sampler.symbol, sampler.context,
               name, name.expression=NULL)}

\arguments{
  \item{sampler.symbol}{A one-element character vector containing
    the symbol of the C sampler function.}
  \item{sampler.context}{An arbitrary R object to be passed to the
    sampler function.}
  \item{name}{The \code{name} attribute for the sampler function.}
  \item{name.expression}{The \code{name.expression} attribute for
    the sampler function, in \code{plotmath} format.}
}

\details{
  This function is intended to allow
  \code{\link{compare.samplers}} to be able to invoke sampler
  functions written in C.  It wraps a C sampler function in an
  R function implementing the standard sampler interface.
  The function named by \code{sampler.symbol} is expected to have
  the function prototype:
  \preformatted{sampler(SEXP sampler_context, dist_t *ds, double *x0,
        int sample_size, double tuning, double *X_out);}
  This is defined as a type \code{sampler_t} in \code{SamplerCompare.h}.
  The parameter \code{sampler_context} is the same as the R object
  \code{sampler.context} passed to \code{wrap.c.sampler}.  \code{ds}
  describes the distribution to be sampled. \code{x0}, \code{sample_size},
  and \code{tuning} are as described in \code{\link{compare.samplers}} and
  should be considered read-only.  \code{X_out} is a column-major
  matrix to be filled in with the generated sample; it has dimension
  \code{ds->ndim * sample_size}.

  The vignette \dQuote{R/C Glue in SamplerCompare} covers this
  interface in greater detail.
}

\value{
  An R function implementing the interface described in
  \code{\link{compare.samplers}}.
}

\seealso{
  \dQuote{R/C Glue in SamplerCompare} (vignette)
}
