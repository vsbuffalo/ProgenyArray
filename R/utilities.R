# utilities.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' Calculate allele frequenceies.
#'
#' @export
alleleFreqs <- function(x) {
  nind <- ncol(x)
  nmissing <- rowSums(is.na(x))
  rowSums(x, na.rm=TRUE)/(2*(nind - nmissing))
}

#' Calculate genotype frequencies.
#'
#' @export
genotypeFreqs <- function(x, counts=FALSE) {
  out <- t(apply(x, 1, countGenotypes))
  if (!counts)
    return(sweep(out, 1, rowSums(tt), "/"))
  return(out)
}

#' Create dataframe of Hardy-Weinberg values.
#'
#' @export
HW <- function(x) {
  out <- as.data.frame(genotypeFreqs(x))
  colnames(out) <- c("g0", "g1", "g2", "gNA")
  out$freq <- alleleFreqs(x)
  out
}

#' Plot genotype and allele frequencies.
#'
#' @export
HWPlot <- function(x) {
  p <- ggplot(HW(x))
  p <- p + geom_point(aes(x=freq, y=g0), color="red")
  p <- p + geom_point(aes(x=freq, y=g1), color="blue")
  p <- p + geom_point(aes(x=freq, y=g2), color="dark green")
  p
}


