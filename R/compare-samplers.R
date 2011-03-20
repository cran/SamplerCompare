# From SamplerCompare, (c) 2010 Madeleine Thompson
# $Id: compare-samplers.R 2756 2011-03-16 11:36:08Z mthompson $

# compare-samplers.R contains compare.samplers and its support
# function eval.sampler, which runs an individual simulation.
# See "Graphical Comparison of MCMC Samplers" (forthcoming) for
# discussion of the figures of merit used here.

# compare.samplers is the main entry point for the SamplerCompare
# package.  See ?compare.samplers for more information.

compare.samplers <- function(sample.size, dists, samplers,
                             tuning=1, trace=TRUE, seed=17,
                             burn.in=0.2, cores=1, completed.file=NULL) {
  # Ensure all distributions are of class dist and all samplers are
  # functions with a name attribute.

  stopifnot(all(sapply(dists, function(a) class(a)=='dist')))

  stopifnot(all(sapply(samplers,
    function(a) class(a)=='function' && !is.null(attr(a,'name')))))

  # Define a temporary save file if necessary.

  if (is.null(completed.file)) {
    unlink.completed.file <- TRUE
    completed.file <- tempfile("compare.samplers")
  } else {
    unlink.completed.file <- FALSE
  }
  if (trace) {
    cat(sprintf('Simulation started at %s.\n', Sys.time()))
    cat(sprintf('Writing results to %s.\n', completed.file))
  }

  # Set up a data frame with a row for each simulation to run.

  jobs <- expand.grid(sampler.id=1:length(samplers), dist.id=1:length(dists),
                      tuning=tuning)

  # Set up locking for an incremental save file.

  if (cores>1) {
    require(synchronicity)
    completed.lock <- boost.mutex()
  } else {
    completed.lock <- NULL
  }
  
  # Save an empty result set with timing attributes to completed.file.

  sampler.comparison <- data.frame()
  attr(sampler.comparison, 'start.time') <- as.numeric(proc.time()['elapsed'])
  attr(sampler.comparison, 'start.stamp') <- Sys.time()
  save(sampler.comparison, file=completed.file)

  # Curry out all the parameters that do not vary simulation to simulation.

  eval.sampler.job.id <- function(job.id) {

    dist <- dists[[jobs$dist.id[job.id]]]
    sampler <- samplers[[jobs$sampler.id[job.id]]]
    tuning.param <- jobs$tuning[job.id]
    eval.sampler(dist, sampler, sample.size, tuning.param,
                 burn.in, trace, seed, completed.file, completed.lock)
  }

  # Call eval.sampler.job.id on each job id, possibly using multiple cores.

  stopifnot(cores>=1)
  if (cores==1) {
    results <- lapply(1:nrow(jobs), eval.sampler.job.id)
  } else {
    require(multicore)
    results <- mclapply(1:nrow(jobs), eval.sampler.job.id,
                        mc.preschedule=FALSE, mc.cores=cores)
  }
#  RS <- do.call(rbind, lapply(results, data.frame))
  
  # Load the saved results and return them.  We could use the lapply
  # return value, but this ensures that callers get the same results
  # whether they get them from compare.samplers or completed.file.

  load(completed.file)
  if (unlink.completed.file)
    unlink(completed.file)
  if (trace)
    cat(sprintf('Simulation finished at %s, %.3gs elapsed.\n',
		 attr(sampler.comparison, 'last.stamp'),
		 attr(sampler.comparison, 'elapsed.time')))
  if (trace && !unlink.completed.file)
    cat(sprintf('Wrote results to %s.\n', completed.file))
	 
  return(sampler.comparison)
}

# Takes a distribution, a sampler, a sample size, and a tuning
# parameter and runs the associated simulation.  Returns a list with
# the distribution name, the sampler name, the tuning parameter, the
# autocorrelation time and information about its uncertainty, the
# number of log density evaluations, the processor time consumed, the
# two-norm of the error in the mean estimate, and a flag indicating
# whether the simulation was aborted.

eval.sampler <- function(dist, sampler, sample.size, tuning.param,
                         burn.in, trace=FALSE, seed=17,
                         completed.file, completed.lock) {

  # Set seed if requested.

  if (seed) {
    if (!exists('.Random.seed'))
      runif(1)
    saved.seed <- .Random.seed
    set.seed(seed)
  }

  # Initialize value to return.

  rv <- list(dist=dist$name, ndim=dist$ndim,
             sampler=attr(sampler, 'name'), tuning=tuning.param,
             act=NA, act.025=NA, act.975=NA,
             act.y=NA, act.y.025=NA, act.y.975=NA,
             evals=NA, grads=NA, cpu=NA, err=NA, aborted=TRUE)

  # Figure out plotmath names for sampler and distribution,
  # defaulting to the regular names.

  a <- attr(sampler, 'name.expression')
  rv$sampler.expr <-
    ifelse(is.null(a), sprintf('plain("%s")', rv$sampler), a)
  a <- dist$name.expression
  rv$dist.expr <-
    ifelse(is.null(a), sprintf('plain("%s")', rv$dist), a)

  # Run the sampler started at a random point on the unit hypercube.

  x0 <- runif(dist$ndim)
  timings <- system.time(
    S <- sampler(target.dist=dist, x0=x0, sample.size=sample.size,
                 tuning=tuning.param) )
  stopifnot(!is.null(S$X) && !is.null(S$evals))

  # Discard the first 20% of the sample and compute the autocorrelation
  # time on what remains.

  discard.obs <- floor(sample.size * burn.in)
  if (nrow(S$X)>discard.obs) {  # enough steps for act
    chain <- S$X[discard.obs:nrow(S$X),,drop=FALSE]
    acts <- ar.act(chain, true.mean=dist$mean)

    rv$act <- acts$act
    rv$act.025 <- acts$act.025
    rv$act.975 <- acts$act.975

    y <- apply(chain, 1, dist$log.density)
    acts.y <- ar.act(y, true.mean=dist$mean.log.dens)
    rv$act.y <- acts.y$act
    rv$act.y.025 <- acts.y$act.025
    rv$act.y.975 <- acts.y$act.975

    # Fill in remaining fields in return value with simulation results.

    rv$evals <- S$evals / nrow(S$X)
    rv$grads <- ifelse(is.null(S$grads), NA, S$grads/nrow(S$X))
    rv$cpu <- timings[['elapsed']] / nrow(S$X)
    rv$aborted <- nrow(S$X) < sample.size

    # Compute the error in the mean if possible.

    if (!is.null(dist$mean))
      rv$err <- sqrt(sum((dist$mean - colMeans(chain))^2))
  }

  # Print results of simulation on terminal if requested.

  if (trace)
    cat(sprintf('%s %s: %.3g (%.3g,%.3g) evals tuning=%.3g%s; act.y=%.3g\n',
                dist$name, attr(sampler,'name'), rv$act * rv$evals, 
                rv$act.025 * rv$evals, rv$act.975 * rv$evals,
                tuning.param, ifelse(rv$aborted, ' (aborted)', ''), rv$act.y))

  # Restore seed if previously modified.

  if (seed)
    .Random.seed <- saved.seed

  # Save results to incremental file if necessary.

  if (!is.null(completed.file)) {
    if (!is.null(completed.lock))
      lock(completed.lock)

    load(completed.file)
    sc.new <- rbind(sampler.comparison, data.frame(rv))
    row.names(sc.new) <- 1:nrow(sc.new)
    attr(sc.new, 'start.time') <- attr(sampler.comparison, 'start.time')
    attr(sc.new, 'start.stamp') <- attr(sampler.comparison, 'start.stamp')
    sampler.comparison <- sc.new
    attr(sampler.comparison, 'elapsed.time') <-
      as.numeric(proc.time()['elapsed']) -
      attr(sampler.comparison, 'start.time')
    attr(sampler.comparison, 'last.stamp') <- Sys.time()
    save(sampler.comparison, file=completed.file)
    if (!is.null(completed.lock))
      unlock(completed.lock)
  }

  return(rv)
}
