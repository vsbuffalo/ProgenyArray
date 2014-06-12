# parentage-methods.R -- statistical methods for inferring parentage
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

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
MendelianTransmissionWithError <- function(ehet=0.6, ehom=0.1) {
  out <- lapply(MendelianTransmissionList(),
                function(x) x %*% genotypingErrorMatrix(ehet, ehom))
  array(unlist(out), dim=c(3, 3, 3))
}

#' Internal function for inferring father using MLE methods
#'
inferFather <- function(progeny, mother, fathers, tmatrix) {
	# first check that fathers is not dimensionless - these means there's only
	# one loci
	if (is.null(dim(fathers)))
		return(list(father=NA, ll=NA, delta=NA, all_lle=NA, nloci=NA))

  # this function gathers all of the probabilities from our transmission matrix,
  # for a particular father, given by their index in fathers.
  tmprob <- function(fi) tmatrix[cbind(fathers[, fi]+1L, progeny+1L, mother+1L)]

  # calculate all of the father's transmission probabilities for each locus
  transmission_probs <- lapply(seq_len(ncol(fathers)), tmprob)

  # calculate log-likelihood under the assumption that all loci are independent
	# FIXME handle NAs
  lle <- sapply(transmission_probs, function(x) sum(log(x)))
  lle_order <- order(lle, decreasing=TRUE)
  delta <- lle[lle_order[1]] - lle[lle_order[2]]
  list(father=lle_order[1], ll=lle[lle_order[1]],
			 delta=delta, all_lle=lle, nloci=length(progeny))
}

#' Check if is valid: all progeny have mothers in genotype matrix.
checkValidProgenyArray <- function(x) {
  if (length(mothers(x)) != ncol(x@progeny_geno))
    stop("mothers vector must be same length as progeny vector")
  if (any(is.na(mothers(x))))
    stop("mothers cannot be NA")
	too_few <- colSums(!is.na(x@parents_geno)) < 100
	if (any(too_few))
		stop(sprintf("%d parents have fewer than 100 loci", sum(too_few))
}


#' Infer Fathers for all progeny when the mother is known
#'
#' @export
setMethod("inferFathers", c(x="ProgenyArray"),
					function(x, ehet, ehom, verbose=FALSE) {
            checkValidProgenyArray(x)
						tmatrix <- MendelianTransmissionWithError(ehet, ehom)
						fathers_lle <- vector('list', ncol(progenyGenotypes((x))))

            # extract the parental and progeny genotypes for complete loci
						parents <- parentGenotypes(x)[x@complete_loci, ]
            kids <- progenyGenotypes(x)[x@complete_loci, ]
            moms_i <- mothers(x)

						for (i in seq_len(ncol(kids))) {
							kid <- kids[, i]
							mom <- parents[, moms_i[i]]
							# find loci that are complete in kids (note all complete in parents)
							noNA_i <- which(!is.na(kid))
							fathers_lle[[i]] <- inferFather(kid[noNA_i], mom[noNA_i],
																							parents[noNA_i, ], tmatrix)
							if (verbose)
								message(sprintf("inferred father of %d (mother = %d) as %d with %d loci",
																i, moms_i[i], fathers_lle[[i]][[1]],
																length(noNA_i)))
						}
						x@fathers <- sapply(fathers_lle, function(x) x[[1]])
            x@fathers_lle <- fathers_lle
						return(x)
})

