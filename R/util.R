# From SamplerCompare, (c) 2010 Madeleine Thompson
# $Id: util.R 1494 2010-08-26 13:30:40Z mthompson $

# util.R is a place where assorted utility functions find a home.

# Not currently used; draws n rows uniformly from the unit sphere in R^p.
# Read http://en.wikipedia.org/wiki/Hypersphere for more.

rsphere <- function(n,p) {
  r <- runif(n)^(1/p)
  x <- array(rnorm(n*p), c(n,p))
  x <- x / sqrt(rowSums(x^2)) * r
  return(x)
}

# Sooner or later, every program I write gets a logsumexp function.
# Sometimes two.

logsumexp <- function(x) {
  m <- max(x)
  return(m+log(sum(exp(x-m))))
}
