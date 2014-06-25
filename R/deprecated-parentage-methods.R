# deprecated-parentage-methods.R -- statistical methods for inferring parentage
# in R (deprecated in favor of the C++/Rcpp methods)
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.
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


##' Infer Fathers for all progeny when the mother is known
##'
#setMethod("inferFathers", c(x="ProgenyArray"),
#          function(x, ehet, ehom, verbose=FALSE) {
#            checkValidProgenyArray(x)
#            fathers_lle <- vector('list', ncol(progenyGenotypes((x))))
#
#            # extract the parental and progeny genotypes for complete loci
#            parents <- parentGenotypes(x)[x@complete_loci, ]
#            kids <- progenyGenotypes(x)[x@complete_loci, ]
#            moms_i <- mothers(x)
#
#            for (i in seq_len(ncol(kids))) {
#              kid <- kids[, i]
#              mom <- parents[, moms_i[i]]
#              # find loci that are complete in kids (note all complete in parents)
#              noNA_i <- which(!is.na(kid))
#              fathers_lle[[i]] <- inferFather(kid[noNA_i], mom[noNA_i],
#                                              parents[noNA_i, ], ehet, ehom)
#              if (verbose)
#                message(sprintf("inferred father of %s (mother = %s) is %s with %d loci",
#                                progenyNames(x)[i], parentNames(x)[moms_i[i]],
#                                parentNames(x)[fathers_lle[[i]][[1]]],
#                                length(noNA_i)))
#            }
#            x@fathers <- sapply(fathers_lle, function(x) x[[1]])
#            #x@lrt <- sapply(fathers_lle, function(x) x[[2]])
#            x@fathers_lle <- fathers_lle
#            names(x@fathers_lle) <- progenyNames(x)
#            return(x)
#})
