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

#' Internal function for inferring father using MLE methods
#'
inferFather <- function(progeny, mother, fathers, ehet, ehom) {
	tmatrix <- MendelianTransmissionWithError(ehet, ehom)
	# first check that fathers is not dimensionless - these means there's only
	# one loci
	if (is.null(dim(fathers)))
		return(list(father=NA, lrt=NA, ll=NA, delta=NA, all_ll=NA, nloci=NA))

  # this function gathers all of the probabilities from our transmission matrix,
  # for a particular father, given by their index in fathers.
  tmprob <- function(fi) tmatrix[cbind(fathers[, fi]+1L, progeny+1L, mother+1L)]

  # calculate all of the father's transmission probabilities for each locus
  transmission_probs <- lapply(seq_len(ncol(fathers)), tmprob)

  # calculate log-likelihood under the assumption that all loci are independent
  ll <- lapply(transmission_probs, function(x) log(x))
	lle <- sapply(ll, sum)
  lle_order <- order(lle, decreasing=TRUE)
	mle_dad <- lle_order[1]
  delta <- lle[mle_dad] - lle[lle_order[2]]

	# psueod-LRT for mle dad
	allele_freqs <- alleleFreqs(cbind(progeny, mother, fathers), min=TRUE)
	stopifnot(min(allele_freqs) > 0)
	# conditional prob of offspring given alleged father
	cpoa <- probOffspringGivenParent(progeny, fathers[mle_dad], allele_freqs, ehet, ehom)
	lo <- log(cpoa)
	lrt <- sum(ll[[mle_dad]] - lo)
  list(father=mle_dad, lrt=lrt, ll=lle[lle_order[1]],
			 delta=delta, all_ll=lle, nloci=length(progeny))
}

inferFatherWithoutMother <- function(progeny, fathers, ehet, ehom) {
	
	cpoa <- probOffspringGivenParent(progeny, fathers, allele_freqs, ehet, ehom)
}

#' Check if is valid: all progeny have mothers in genotype matrix.
checkValidProgenyArray <- function(x) {
  if (length(mothers(x)) != ncol(x@progeny_geno))
    stop("mothers vector must be same length as progeny vector")
  if (any(is.na(mothers(x))))
    stop("mothers cannot be NA")
	too_few <- colSums(!is.na(x@parents_geno)) < 100
	if (any(too_few))
		stop(sprintf("%d parents have fewer than 100 loci", sum(too_few)))
}


#' Infer Fathers for all progeny when the mother is known
#'
#' @export
setMethod("inferFathers", c(x="ProgenyArray"),
					function(x, ehet, ehom, verbose=FALSE) {
            checkValidProgenyArray(x)
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
																							parents[noNA_i, ], ehet, ehom)
							if (verbose)
								message(sprintf("inferred father of %s (mother = %s) is %s with %d loci",
																progenyNames(x)[i], parentNames(x)[moms_i[i]],
															 	parentNames(x)[fathers_lle[[i]][[1]]],
																length(noNA_i)))
						}
						x@fathers <- sapply(fathers_lle, function(x) x[[1]])
						#x@lrt <- sapply(fathers_lle, function(x) x[[2]])
            x@fathers_lle <- fathers_lle
						names(x@fathers_lle) <- progenyNames(x)
						return(x)
})

setMethod("parentage", c(x="ProgenyArray"),
					function(x) {
						data.frame(progeny=progenyNames(x),
											 mother=parentNames(x)[mothers(x)],
											 father=parentNames(x)[fathers(x)])

					})
