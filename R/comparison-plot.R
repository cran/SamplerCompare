# From SamplerCompare, (c) 2010 Madeleine Thompson
# $Id: comparison-plot.R 3045 2011-04-23 03:03:14Z mthompson $

# comparison.plot generates a plot comparing MCMC performance as
# described in "Graphical Comparison of MCMC Performance" (forthcoming).
# See ?comparison.plot for usage details.

comparison.plot <- function(RS, xlab=NULL, ylab=NULL, base_size=10, ...) {
  stopifnot(require(ggplot2))

  # First, plot the results with finite ACT.

  RSfinite <- subset(RS, is.finite(RS$act))

  # Compute a reasonably-spaced set of tick marks.

  x.breaks <- log.breaks(RSfinite$tuning, 10)

  # Manually handle xlab and ylab so they can be overridden by caller.

  if (is.null(xlab))
    xlab <- 'scale tuning parameter'
  if (is.null(ylab))
    ylab <- paste('# of evals. of log density function per',
                  'uncorrelated obs. (with 95% CI)')

  # Generate grid of plots of evals*act vs. tuning parameter, with
  # a 95% confidence interval.  The theme and x scale overrides are
  # largely personal preference; if the caller wants different choices,
  # they can specify them with the + operator on the returned value.
  # "labeller" really is spelled this way in facet_grid.

  p <- qplot(tuning, evals*act, ymin=evals*act.025, ymax=evals*act.975,
      data=RSfinite, log='xy', geom='pointrange', xlab=xlab, ylab=ylab, ...) +
    facet_grid(dist.expr~sampler.expr, labeller=label_parsed) +
    scale_x_log(breaks=x.breaks, labels=x.breaks) +
    theme_bw(base_size=base_size) +
    opts(panel.grid.minor=theme_blank(),
         axis.text.x=theme_text(angle=45, vjust=1))

  # Next, plot the results with infinite/unknown ACT as question
  # marks at lim, computed to be at the top of the plot.
  
  if (any(!is.finite(RS$act))) {
    RSinf <- subset(RS, !is.finite(RS$act))
    RSinf$lim <- max(RSfinite$evals*RSfinite$act)
    p <- p + geom_text(aes(x=tuning, y=lim, label='?'), data=RSinf, size=4.0)
  }

  return(p)
}

# Return an sequence of tick marks to be used on a log axis.  Used
# to keep the automatic tick generation from spacing ticks too
# closely; this only returns ticks on integer powers of base.

log.breaks <- function(data, base) {
  breaks <- base^seq(floor(min(log(data,base=base))),
                     ceiling(max(log(data,base=base))))
}
