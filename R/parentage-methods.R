# parentage-methods.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' which.max for matrices
whichRowColMax <- function(x) {
  c(row(x)[which.max(x)], col(x)[which.max(x)])
}

#' calculate LODs for QQ/UU and QU/UU, for maximum parent 1 and 2
caclulateLODs <- function(x) {
  mlparents <- lapply(x, function(l) whichRowColMax(l[[1]] - l[[2]]))
  #mlparents <- lapply(x, function(l) whichRowColMax(l[[1]]))
  # using Meagher and Thompson 1985 comparisons
  qu_comparisons <- do.call(rbind, lapply(seq_along(mlparents), function(i) {
                            pars <- mlparents[[i]]
                            qq <- x[[i]]$prob_qq
                            uu <- x[[i]]$prob_uu
                            data.frame(qq_uu=(qq - uu)[pars[1], pars[2]])
                           }))
  d <- qu_comparisons
  d
}

#' Infer parents of all progeny in a ProgenyArray
#'
#' inferParents() infers the parents of all ProgenyArray by calculating the
#' likelihood of all possible parents. If there are P parents, this involves
#' calculating P(P - 1)/2 likelihoods.
#' @param x a ProgenyArray object
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#' @param verbose logical indicating whether to use verbose output
#' @export
setMethod("inferParents", c(x="ProgenyArray"),
function(x, ehet, ehom, verbose=TRUE) {
  #freqs <- alleleFreqs(progenyGenotypes(x))[x@complete_loci]
  freqs <- freqs(pa)[x@complete_loci]
  stopifnot(all(freqs > 0 & freqs < 1)) # nothing fixed
  parents <- parentGenotypes(x)[x@complete_loci, ]
  kids <- progenyGenotypes(x)[x@complete_loci, ]
  pars <- vector("list", ncol(kids)) # for ths data from .allParentLikelihoods
  out <- list()

  if (verbose)
    message(sprintf("inferring parents for %d progeny", ncol(kids)))

  pars <- .inferParents(kids, parents, freqs, ehet, ehom, verbose)

  # calculate the LOD scores for parents, extract the ML parent
  lods <- caclulateLODs(pars)
  mlparents <- lapply(pars, function(l) whichRowColMax(l[[1]] - l[[2]]))
  #mlparents <- lapply(pars, function(l) whichRowColMax(l[[1]]))

  nloci <- sapply(pars, '[[', 3)

  # use parents slot, since we don't know who mom and dad are
  x@parents <- data.frame(progeny=seq_len(ncol(kids)),
                          parent_1=sapply(mlparents, '[', 1),
                          parent_2=sapply(mlparents, '[', 2),
                          lods=lods)
 # find given mothers that are inconsistent with the one found
  if (length(x@supplied_mothers)) {
    mothers <- x@supplied_mothers
    moms <- sapply(seq_along(mlparents), function(i) {
                   match(mothers[i], mlparents[[i]])
                  })
    inconsistent_moms <- which(is.na(moms))
    #stopifnot(length(moms) == length(x@supplied_mothers))
    ninc <- length(inconsistent_moms)
    if (ninc > 0L)
      warning(sprintf("found %d mothers that are inconsistent", ninc))
    x@parents$which_mother <- moms

  }

  # add debugging information
  debug <- data.frame(nloci=nloci)
  x@parents <- cbind(x@parents, debug)
  x@parent_lods <- pars
  return(x)
})

