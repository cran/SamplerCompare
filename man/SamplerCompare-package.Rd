% From SamplerCompare, (c) 2010-2011 Madeleine Thompson

\name{SamplerCompare-package}
\alias{SamplerCompare-package}
\alias{SamplerCompare}
\docType{package}
\title{A framework for comparing the performance of MCMC samplers}
\description{This package consists of two components: a framework
for running sets of MCMC samplers on sets of distributions with a
variety of tuning parameters and plotting functions to visualize
the results of those simulations.  See sc-intro.pdf for an introduction.}
\details{
\tabular{ll}{
Package: \tab SamplerCompare\cr
Type: \tab Package\cr
Version: \tab 1.2.4\cr
Date: \tab 2013-04-18\cr
License: \tab GPL-2\cr
Copyright: \tab (c) 2010-2013 Madeleine Thompson, except \file{dchud.f} and
  \file{dchdd.f}, which are in the public domain in the United States.\cr
LazyLoad: \tab yes\cr
}

The central function of this package is \code{\link{compare.samplers}},
which simulates a set of distributions with a set of samplers with
a set of tuning parameters.  It returns a data frame which may be
passed to \code{\link{comparison.plot}}.

The source files \file{src/dchud.f} and \file{src/dchdd.f} come from LINPACK
and were written by G. W. Stewart.

The rest of the package was written by Madeleine Thompson.
}
\author{Madeleine Thompson, except dchud.f and dchdd.f, which were written by
  G. W. Stewart.}
