# parentage-methods.R -- statistical methods for inferring parentage
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#MT <- matrix(c(0, 0.5, 1), ncol=1)

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

# convert to array. Now dimensions are:
# (1) father, (2) progeny, (3) mother
MT <- array(unlist(MT_list), dim=c(3, 3, 3))
stopifnot(all(sapply(1:3, function(i) MT_list[[i]] == MT[, , i])))

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
  matrix(c(1-ehom, ehet/2, ehom/2, ehom/2, 1-ehom, ehom/2, ehom/2, ehet/2, 1-ehom),
         ncol=3)
}


#' Create a Mendelian transmission matrix, with error in observed progeny genotypes.
#'
#' @param ehet the probability that a heterozygous genotype is incorrect
#' @param ehom the probability that a homozygous genotype is incorrect
#' @export
transmissionMatrix <- function(ehet=0.6, ehom=0.1) {
  out <- lapply(MT_list, function(x) x %*% genotypingErrorMatrix(ehet, ehom))
  array(unlist(out), dim=c(3, 3, 3))
}

#' Internal function for inferring father using MLE methods
#'
inferFather <- function(progeny, mother, fathers, tmatrix) {
  # this function gathers all of the probabilities from our transmission matrix,
  # for a particular father, given by their index in fathers.
  tmprob <- function(fi) tmatrix[cbind(fathers[, fi]+1L, progeny+1L, mother+1L)]

  # calculate all of the father's transmission probabilities for each locus
  transmission_probs <- lapply(seq_len(ncol(fathers)), tmprob)

  # calculate log-likelihood under the assumption that all loci are independent
  lle <- sapply(transmission_probs, function(x) sum(log(x)))
  lle_order <- order(lle, decreasing=TRUE)
  delta <- lle[lle_order[1]] - lle[lle_order[2]]
  list(father=lle_order[1], ll=lle[lle_order[1]], delta=delta, all_lle=lle)
}

#' Infer Fathers for all progeny when the mother is known
#'
#' @export
setMethod("inferFathers", c(x="ProgenyArray"),
					function(x, ehet, ehom) {
						tmatrix <- transmissionMatrix(ehet, ehom)
})

