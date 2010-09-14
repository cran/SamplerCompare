# From SamplerCompare, (c) 2010 Madeleine Thompson
# $Id: compare-samplers.R 1603 2010-09-11 13:51:42Z mthompson $

# compare-samplers.R contains compare.samplers and its support
# function eval.sampler, which runs an individual simulation.
# See "Graphical Comparison of MCMC Samplers" (forthcoming) for
# discussion of the figures of merit used here.

# compare.samplers is the main entry point for the SamplerCompare
# package.  See ?compare.samplers for more information.

compare.samplers <- function(sample.size, dists, samplers,
                             tuning=1, trace=TRUE, seed=17,
                             include.ics=FALSE) {
  # Ensure all distributions are of class dist and all samplers are
  # functions with a name attribute.

  stopifnot(all(sapply(dists, function(a) class(a)=='dist')))

  stopifnot(all(sapply(samplers,
    function(a) class(a)=='function' && !is.null(attr(a,'name')))))

  results <- list()

  # Iterate over distributions, samplers, and tuning parameters.

  for (dist in dists) {
    for (sampler in samplers) {
      for (tuning.param in tuning) {

        # Set seed if requested.

        if (seed) {
          if (!exists('.Random.seed'))
            runif(1)
          saved.seed <- .Random.seed
          set.seed(seed)
        }

        # Run one simulation.
        
        run <- eval.sampler(dist, sampler, sample.size, tuning.param,
                            include.ics)

        # Print results of simulation on terminal if requested.

        if (trace)
          cat(sprintf('%s %s: %.3g (%.3g,%.3g) evals tuning=%.3g %s\n',
                      dist$name, attr(sampler,'name'), run$act * run$evals, 
                      run$act.025 * run$evals, run$act.975 * run$evals,
                      tuning.param, ifelse(run$aborted, '(aborted)', '')))

	# Figure out plotmath names for sampler and distribution,
	# defaulting to the regular names.

        a <- attr(sampler, 'name.expression')
        run$sampler.expr <-
          ifelse(is.null(a), sprintf('plain("%s")', run$sampler), a)
        a <- dist$name.expression
        run$dist.expr <-
          ifelse(is.null(a), sprintf('plain("%s")', run$dist), a)

        # Add results from this run to the result set.

        results <- rbind(results, run)

        # Restore seed if previously modified.
        if (seed)
          .Random.seed <- saved.seed
      }
    }
  }
  
  # This is not obvious.  First, turn the results list into a data
  # frame with lists for each column.  Use lapply with unlist to turn
  # it into a list of primitive vectors/factors.  Then use data.frame
  # to turn it back into a data frame with vectors for each column.
  # The row names must be set to NULL or else the data.frame calls
  # will attempt to infer row names, and decide that since each row
  # came from a same-named variable, all the row names are the same,
  # which results in a spurious warning.  It seems like there should
  # be a simpler way to do this.

  row.names(results) <- NULL
  RS <- data.frame(lapply(data.frame(results),unlist))
  attr(RS, 'sample.size') <- sample.size
  return(RS)
}

# Takes a distribution, a sampler, a sample size, and a tuning
# parameter and runs the associated simulation.  Returns a list with
# the distribution name, the sampler name, the tuning parameter, the
# autocorrelation time and information about its uncertainty, the
# number of log density evaluations, the processor time consumed, the
# two-norm of the error in the mean estimate, and a flag indicating
# whether the simulation was aborted.

eval.sampler <- function(dist, sampler, sample.size, tuning.param,
                         include.ics) {

  # Initialize value to return.

  rv <- list(dist=dist$name, ndim=dist$ndim,
             sampler=attr(sampler, 'name'), tuning=tuning.param,
             act=NA, act.025=NA, act.975=NA, act.se=NA,
             evals=NA, grads=NA, cpu=NA, err=NA, aborted=TRUE)

  # Run the sampler started at a random point on the unit hypercube.

  x0 <- runif(dist$ndim)
  timings <- system.time(
    S <- sampler(target.dist=dist, x0=x0, sample.size=sample.size,
                 tuning=tuning.param) )
  stopifnot(!is.null(S$X) && !is.null(S$evals))

  # Discard the first 20% of the sample and compute the autocorrelation
  # time on what remains.

  burnin <- floor(sample.size * 0.2)
  if (nrow(S$X)<=burnin)  # not enough steps for act
    return(rv)
  chain <- S$X[burnin:nrow(S$X),,drop=FALSE]
  acts <- ar.act(chain, true.mean=dist$mean)

  rv$act <- acts$act
  rv$act.025 <- acts$act.025
  rv$act.975 <- acts$act.975
  rv$act.se <- acts$se

  # If requested, use the initial convex sequence estimator of C. Geyer (1992).

  if (include.ics)
    rv$act.ics <- ics.act(chain, true.mean=dist$mean)

  # Fill in remaining fields in return value with simulation results.

  rv$evals <- S$evals / nrow(S$X)
  rv$grads <- ifelse(is.null(S$grads), NA, S$grads/nrow(S$X))
  rv$cpu <- timings[['elapsed']] / nrow(S$X)
  rv$aborted <- nrow(S$X) < sample.size

  # Compute the error in the mean if possible.

  if (!is.null(dist$mean))
    rv$err <- sqrt(sum((dist$mean - colMeans(chain))^2))

  return(rv)
}
