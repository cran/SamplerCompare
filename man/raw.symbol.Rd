% From SamplerCompare, (c) 2010 Madeleine Thompson
% $Id: raw.symbol.Rd 1487 2010-08-25 21:48:40Z mthompson $

\name{raw.symbol}

\alias{raw.symbol}

\title{Locate a symbol}

\description{Call R_FindSymbol and return function pointer in a raw vector}

\usage{raw.symbol(symbol)}

\arguments{
  \item{symbol}{a length one character vector containing a C symbol}
}

\details{
  This function calls \code{R_FindSymbol(symbol, "", NULL)} in C.
  If the symbol is found, the function pointer is returned as a raw
  vector.  If not, an error is thrown.

  This is intended to be used to fill in context objects for samplers
  and distributions implemented in C.  Exposing this interface in
  R prevents the need for extra C glue that does nothing except
  call \code{R_FindSymbol}.
}

\value{
  A raw vector containing a function pointer.
}

\seealso{
  \code{\link{wrap.c.sampler}},
  \code{\link{make.c.dist}},
  \dQuote{R/C Glue in SamplerCompare} (vignette)
}
