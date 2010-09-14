# From SamplerCompare, (c) 2010 Madeleine Thompson
# $Id: dist-util.R 1494 2010-08-26 13:30:40Z mthompson $

# This file contains code for managing objects of the "dist" class,
# which represent probability distributions in the SamplerCompare
# package.  The R help for these functions and the vignette "R/C Glue
# in SamplerCompare (doc/glue.pdf) may also be of interest.

# See ?make.dist for more information.

make.dist <- function(ndim, name, name.expression=NULL,
                      log.density=NULL, grad.log.density=NULL,
                      log.density.and.grad=NULL, mean=NULL, cov=NULL) {
  stopifnot(ndim>0)
  if (!is.null(mean))
    stopifnot(length(mean)==ndim)
  if (!is.null(cov))
    stopifnot(all(dim(cov)==ndim))

  # Define a dist object as a list, filling in whatever the user gave us.

  ds <- list(ndim=ndim, name=name, name.expression=name.expression,
             log.density=log.density, grad.log.density=grad.log.density,
             log.density.and.grad=log.density.and.grad, mean=mean, cov=cov)

  # If the user gave us log.density.and.grad but not log.density,
  # fake up a log.density.

  if (is.null(log.density) && !is.null(log.density.and.grad)) {
    ds$log.density <- function(x) log.density.and.grad(x,FALSE)$log.density
  }

  # If the user gave us log.density.and.grad but not grad.log.density,
  # fake up a grad.log.density.

  if (is.null(grad.log.density) && !is.null(log.density.and.grad)) {
    ds$grad.log.density <-
      function(x) log.density.and.grad(x,TRUE)$grad.log.density
  }

  # If the user gave us log.density and grad.log.density but not
  # log.density.and.grad, fake up a log.density.and.grad.

  if (is.null(log.density.and.grad) && !is.null(log.density) &&
      !is.null(grad.log.density)) {
    ds$log.density.and.grad <- function(x, compute.grad=FALSE) {
      if (compute.grad) {
        return(list(log.density=ds$log.density(x),
                    grad.log.density=ds$grad.log.density(x)))
      } else {
        return(list(log.density=ds$log.density(x)))
      }
    }
  }

  # Mark the distribution with its class and return it.

  class(ds) <- 'dist'
  return(ds)
}

# Overrides print() for dist objects.  Prints basic information
# about a distribution.

print.dist <- function(x, ...) {
  if (!is.null(x$c.log.density.and.grad))
    feat <- 'log-density and gradient implemented in C'
  else if (!is.null(x$log.density) && !is.null(x$grad.log.density))
    feat <- 'known log-density and gradient'
  else if (!is.null(x$log.density))
    feat <- 'known log-density'
  else
    feat <- 'unknown log-density'
  cat(sprintf("%s (%d dimensions, %s)\n", x$name, x$ndim, feat))
}

# See ?make.c.dist for more information.

make.c.dist <- function(ndim, name, c.log.density, c.context=NULL,
                        name.expression=NULL, mean=NULL, cov=NULL) {

  # Make a distribution object.

  ds <- make.dist(ndim, name, name.expression=name.expression,
                  mean=mean, cov=cov)
  
  # Locate the symbol for the log density.

  ndim <- as.integer(ndim)
  stopifnot(is.character(c.log.density) && length(c.log.density==1))
  ds$sym <- raw.symbol(c.log.density)

  # Define log.density.and.grad wrapper that calls the C version.

  ds$log.density.and.grad <- function(x, compute.grad=FALSE) {
    stopifnot(is.numeric(x) && length(x)==ndim)
    r <- .Call(R_invoked_C_glue, ds$sym, c.context, x, compute.grad)
    return(r)
  }

  # Define log.density and grad.log.density as wrappers.

  ds$log.density <- function(x) {
    ds$log.density.and.grad(x)$log.density
  }

  ds$grad.log.density <- function(x) {
    ds$log.density.and.grad(x, compute.grad=TRUE)$grad.log.density
  }

  # Note c.log.density.and.grad and c.context for debugging and print().

  ds$c.log.density.and.grad <- c.log.density
  ds$c.context <- c.context

  return(ds)
}

# Ensures that the log density and its gradient in a distribution
# object are consistent with each other at point x.  Uses the central
# difference formula on each coordinate in turn.
#
# See ?check.dist.gradient for more information.  ?make.dist also has
# an example.  See Nocedal & Wright, Numerical Optimization 2nd ed.
# sec. 8.1 (p. 196-197) for more information on the central difference
# formula.

check.dist.gradient <- function(ds, x, h=1e-7) {

  # Compute gradient with analytic formula.

  g <- ds$grad.log.density(x)
  stopifnot(all(is.finite(g)))
  stopifnot(length(g)==length(x))

  for (i in 1:length(x)) {

    # Check analytic gradient against numerical derivative in dimension i.

    e <- c(rep(0,i-1),1,rep(0,length(x)-i))               # element vector
    g.discrete <- (ds$log.density(x+e*h*x[i]) - ds$log.density(x-e*h*x[i])) /
      2 / h / x[i]  # approximation
    stopifnot(all(is.finite(g.discrete)))

    # If relative error is more than 0.001, fail.

    if(abs( (g.discrete-g[i])/g[i] ) > 1e-3) {
      stop("Gradients do not match.\n",
           sprintf("g[%d] = %.5g   g.discrete[%d] = %.5g\n",
                   i, g[i], i, g.discrete))
    }
  }
}
