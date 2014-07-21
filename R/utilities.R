# utilities.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

validGenotypes <- function(x, allowNA=TRUE) {
  genos <- seq(0, 2)
  if (allowNA)
    genos <- c(NA, genos)
  all(x %in% genos)
}


#' Create a genotyping error matrix.
#'
#' Returns a matrix of transition probabilities for a given genotype given error.
#'
#' @param ehet the probability that a heterozygous genotype is incorrect
#' @param ehom the probability that a homozygous genotype is incorrect
#'
#' This function creates a matrix of the form:
#'   [1-e, e/2, e/2,
#'    E/2, 1-E, E/2,
#'    e/2, e/2, 1-e]
#' where E is the probability the heterozygous genotype call is incorrect, and
#' e is the probability that the homozygous genotypes call is incorrect.
#' @export
genotypingErrorMatrix <- function(ehet=0.6, ehom=0.1) {
  matrix(c(1-ehom, ehet/2, ehom/2, ehom/2, 1-ehet, ehom/2, ehom/2, ehet/2, 1-ehom),
         ncol=3)
}

#' Calculate allele frequenceies.
#'
#' @param x a genotype matrix
#' @param min force minimum allele frequency to be 1/nsample
#' @export
alleleFreqs <- function(x, min=FALSE) {
  # TODO make this handle missing data better.
  nind <- ncol(x)
  nmissing <- rowSums(is.na(x))
  out <- rowSums(x, na.rm=TRUE)/(2*(nind - nmissing))
  if (min)
    return(ifelse(out == 0, 1/ncol(x), out))
  return(out)
}

#' Calculate genotype frequencies.
#'
#' @param x a genotype matrix
#' @param counts whether to return genotype counts rather than frequency
#' @export
genotypeFreqs <- function(x, counts=FALSE) {
  out <- t(apply(x, 1, countGenotypes))
  if (!counts)
    return(sweep(out, 1, rowSums(tt), "/"))
  return(out)
}

#' Create dataframe of Hardy-Weinberg values.
#'
#' @param x a genotype matrix
#' @export
HW <- function(x) {
  out <- as.data.frame(genotypeFreqs(x))
  colnames(out) <- c("g0", "g1", "g2", "gNA")
  out$freq <- alleleFreqs(x)
  out
}

#' Plot genotype and allele frequencies.
#'
#' @param x a genotype matrix
#' @export
HWPlot <- function(x) {
  p <- ggplot(HW(x))
  p <- p + geom_point(aes(x=freq, y=g0), color="red")
  p <- p + geom_point(aes(x=freq, y=g1), color="blue")
  p <- p + geom_point(aes(x=freq, y=g2), color="dark green")
  p
}


MendelianTransmissionList <- function(x) {
    # Mendelian transition matrix, in list notation for readability
    # TODO: we can build this up programmatically using Kronecker products.
    # first dimension (list element): mom's genotype
    # second dimension: alleged dad's genotype
    # third dimension: progeny genotype
    MT_list <- list(
               "0"=matrix(c(1, 1/2, 0, 0, 1/2, 1, 0, 0, 0),
                          ncol=3),
               "1"=matrix(c(1/2, 1/4, 0, 1/2, 1/2, 1/2, 0, 1/4, 1/2),
                          ncol=3),
               "2"=matrix(c(0, 0, 0, 1, 1/2, 0, 0, 1/2, 1),
                          ncol=3))
    return(MT_list)
}

MendelianTransmissionMatrix <- function() {
  # convert to array. Now dimensions are:
  # (1) father, (2) progeny, (3) mother
  m <- MendelianTransmissionList()
  MT <- array(unlist(m), dim=c(3, 3, 3))
  stopifnot(all(sapply(1:3, function(i) m[[i]] == MT[, , i])))
  return(MT)
}

conditionalOffspringParentMatrix <- function(freq) {
  #  P(Go | Gm)
  #      offspring
  #    0   1          2
  # 0  b   c          0
  # 1  b/2  (b+c)/2   c/2
  # 2  0    b         c
  b <- 1 - freq
  c <- freq
  matrix(c(b, b/2, 0, c, (b+c)/2, b, 0, c/2, c), nrow=3)
}


#' Create a Mendelian transmission matrix, with error in observed progeny genotypes.
#'
#' @param ehet the probability that a heterozygous genotype is incorrect
#' @param ehom the probability that a homozygous genotype is incorrect
#' @export
MendelianTransmissionWithError <- function(ehet=0.6, ehom=0.1) {
  out <- lapply(MendelianTransmissionList(),
                function(x) x %*% genotypingErrorMatrix(ehet, ehom))
  array(unlist(out), dim=c(3, 3, 3))
}

probOffspringGivenParent <- function(offspring, parent, freqs, ehet, ehom) {
  stopifnot(!any(is.na(offspring)))
  stopifnot(!any(is.na(parent)))
  stopifnot(!any(is.na(freqs)))
  error_matrix <- genotypingErrorMatrix(ehet, ehom)
  mapply(function(f, o, m) {
         mat <- conditionalOffspringParentMatrix(f)
         (mat %*% error_matrix)[m+1L, o+1L]
  }, freqs, offspring, parent)
}


MendelianInconsistencies <- function(parent1, child, parent2=NULL) {
  m <- MendelianTransmissionMatrix() > 0
  stopifnot(validGenotypes(c(parent1, child)))
  parent1 <- parent1+1L
  child <- child+1L
  if (is.null(parent2)) {
    validGenotypes(parent2)
    return(mapply(function(p1, c) all(!m[, c, p1]), parent1, child))
  }
  parent2 <- parent2+1L
  return(mapply(function(p1, p2, c) !m[p2, c, p1], parent1, parent2, child))
}
