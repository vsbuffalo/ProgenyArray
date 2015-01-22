# progenyarray-methods.R -- method for working with ProgenyArray objects
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' Return which loci are complete
whichLociComplete <- function(x) {
  which(apply(!is.na(x), 1, all))
}

#' Which loci are not fixed
whichLociPolymorphic <- function(x) {
  freq <- alleleFreqs(x)
  which(freq > 0 & freq < 1)
}

#' Manually set the allele frequencies
#'
#' @param x a ProgenyArray object
#'
#' @name progenyNames
#' @export
setReplaceMethod("freqs", "ProgenyArray", function(x, value) {
  x@freqs <- value
  # prune fixed or incomplete loci
  not_fixed <- whichLociPolymorphic(x@parents_geno)
  x@complete_loci <- intersect(whichLociComplete(x@progeny_geno), not_fixed)
  x
})


#' Constructor new ProgenyArray object.
#'
#' Create a new ProgenyArray object for storing genotyping and other data from
#' a progeny array experiment.
#'
#' @param progeny_geno progeny genotypes matrix (all entries should be 0, 1, 2, or NA)
#' @param parents_geno parents genotypes matrix (all entries should be 0, 1, 2, or NA)
#' @param mothers user-supplied mothers (e.g. known from collection)
#' @param loci a \code{GRanges} object of loci
#' @param ref a character vector of reference loci
#' @param alt a character vector of alternate loci
#' @export
ProgenyArray <- function(progeny_geno, parents_geno, mothers=integer(),
                         loci=GRanges(), ref=character(), alt=character()) {
  obj <- new("ProgenyArray", progeny_geno=progeny_geno, parents_geno=parents_geno,
             supplied_mothers=mothers, ranges=loci, ref=ref, alt=alt)

  obj@complete_loci <- whichLociComplete(obj@parents_geno)
  not_fixed <- whichLociPolymorphic(obj@parents_geno)
  message(sprintf("%d loci are fixed", nrow(obj@progeny_geno) - length(not_fixed)))
  obj@complete_loci <- intersect(obj@complete_loci, not_fixed)
  obj
}


#' Pretty-print a ProgenyArray object
#'
#' @param object a ProgenyArray object
#'
#' @export
setMethod("show",
          c(object="ProgenyArray"),
          function(object) {
            cat(sprintf("ProgenyArray object: %d loci, %d parents, %d progeny\n",
                        nrow(object@progeny_geno), ncol(object@parents_geno),
                        ncol(object@progeny_geno)))
            cat(sprintf(" number of chromosomes: %d\n object size: %s Mb\n",
                        length(seqlevels(object@ranges)),
                        round(object.size(object)/1024^2, 3)))
            # convience function for getting lengths
            numOrNA <- function(x, unique=FALSE) {
              if (length(x) == 0)
                return(NA)
              if (is.null(dim(x))) {
                if (unique) return(length(unique(x)))
                return(length(x))
              }
              return(ncol(x))
            }
            #nfathers <- numOrNA(object@fathers)
            nparents <- numOrNA(object@parents_geno)
            #nmothers <- numOrNA(object@mothers, unique=TRUE)
            nprogeny <- numOrNA(object@progeny_geno)
            cat(sprintf(" number of progeny: %d\n", nprogeny))
            cat(sprintf(" number of parents: %d\n", nparents))
            #cat(sprintf("Number of fathers: %d\n", nfathers))
            #cat(sprintf("Number of mothers: %d\n", nmothers))
            cat(sprintf(" proportion missing:\n   progeny: %0.3f\n  parents:  %0.3f\n",
                        sum(is.na(object@progeny_geno))/length(object@progeny_geno),
                        sum(is.na(object@parents_geno))/length(object@parents_geno)))
            cat(sprintf(" number of complete parental loci: %d\n",
                        length(object@complete_loci)))

          })


#' Accessor for parent genotypes
#'
#' @param x a ProgenyArray object
#' @export
setMethod("parentGenotypes",
          c(x="ProgenyArray"),
          function(x, seqname=NULL) {
            if (is.null(seqname))
              return(x@parents_geno)
            else
              return(x@parents_geno[as.logical(seqnames(x@ranges) == seqname), ])
          })

#' Accessor for progeny genotypes
#'
#' @param x a ProgenyArray object
#' @export
setMethod("progenyGenotypes",
          c(x="ProgenyArray"),
          function(x, seqname=NULL) {
            if (is.null(seqname))
              return(x@progeny_geno)
            else
              return(x@progeny_geno[as.logical(seqnames(x@ranges) == seqname), ])
          })

#' Accessor for parents \code{data.frame} in a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("parents",
          c(x="ProgenyArray"),
          function(x) {
            return(x@parents)
          })

#' Method to check if individual is selfed
#'
#' @param x a ProgenyArray object
#' @export
setMethod("isSelfed",
          c(x="ProgenyArray"),
          function(x) {
            if (!length(x@supplied_mothers))
              stop("this method only works with user-supplied mothers")
            mothers(x) == fathers(x)
          })


#' Accessor for user-supplied mothers
#'
#' @param x a ProgenyArray object
#' @export
setMethod("suppliedMothers",
          c(x="ProgenyArray"),
          function(x) {
            return(x@supplied_mothers)
          })

#' Accessor for fathers in a ProgenyArray object
#'
#' Internally, this references the \code{parents} slot.
#' @param x a ProgenyArray object
#' @export
setMethod("fathers",
          c(x="ProgenyArray"),
          function(x) {
            dd <- parents(x)
            ifelse(dd$which_mother == 2, dd$parent_1, dd$parent_2)
          })

#' Accessor for mothers in a ProgenyArray object
#'
#' Internally, this references the \code{parents} slot.
#' @param x a ProgenyArray object
#' @export
setMethod("mothers",
          c(x="ProgenyArray"),
          function(x) {
            dd <- parents(x)
            ifelse(dd$which_mother == 1, dd$parent_1, dd$parent_2)
          })

#' Accessor for genomic ranges from a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("granges",
          c(x="ProgenyArray"),
          function(x, use.mcols=FALSE, ...) {
            return(x@ranges)
          })

#' Accessor for alterate alleles from a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("alt",
          c(x="ProgenyArray"),
          function(x) {
            return(x@alt)
          })

#' Accessor for reference alleles from a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("ref",
          c(x="ProgenyArray"),
          function(x) {
            return(x@ref)
          })

#' Return progeny sample names
#'
#' @param x a ProgenyArray object
#' @export
setMethod("progenyNames", "ProgenyArray", function(x) {
  colnames(x@progeny_geno)
})

#' Return parent sample names
#'
#' @param x a ProgenyArray object
#' @export
setMethod("parentNames", "ProgenyArray", function(x) {
  colnames(x@parents_geno)
})

#' Set method for supplied mothers
#'
#' @param x a ProgenyArray object
#' @param value vector of mother indices
#' @export
#' @name suppliedMothers
#'
setReplaceMethod("suppliedMothers", "ProgenyArray", function(x, value) {
  if (length(value) != ncol(x@progeny_geno))
    stop("length of value must be same as number of progeny")
  x@supplied_mothers <- value
  return(x)
})


#' Set method for parents sample names
#'
#' @param x a ProgenyArray object
#' @param value vector of parent names
#' @export
#' @name parentNames
setReplaceMethod("parentNames", "ProgenyArray", function(x, value) {
  if (length(value) != ncol(x@parents_geno))
    stop("length of value must be same as number of parents")
  if (any(duplicated(value)))
    stop("no duplicated names allowed")
  colnames(x@parents_geno) <- value
  return(x)
})


#' Set method for progeny sample names
#'
#' @param x a ProgenyArray object
#' @param value vector of parent names
#' @name progenyNames
#' @export
setReplaceMethod("progenyNames", "ProgenyArray", function(x, value) {
  if (length(value) != ncol(x@progeny_geno))
    stop("length of value must be same as number of progeny")
  if (any(duplicated(value)))
    stop("no duplicated names allowed")
  colnames(x@progeny_geno) <- value
  return(x)
})

#' Get inferred parents as a dataframe
#'
#' @param x a ProgenyArray object
#' @name progenyNames
#' @export

setMethod("parentage", "ProgenyArray", function(x) {
  if (!length(x@supplied_mothers)) stop("if user-supplied mothers are set, use parents() method")
  nprogeny <- ncol(x@progeny_geno)
  nparents <- ncol(x@parents_geno)
  data.frame(progeny=progenyNames(x)[seq_len(nprogeny)],
             mothers=parentNames(x)[mothers(x)],
             fathers=parentNames(x)[fathers(x)])
})


#' Create a diagnostic plot of a parentage inference, looking at all log likelihood ratios.
#'
#' @param x a ProgenyArray object
#' @param progeny which progeny to use (an index to its column in progeny genotypes)
#'
#' @export
setMethod("parentLLPlot", "ProgenyArray", function(x, progeny) {
  lods <- x@parent_lods
  d <- melt(lods[[progeny]][[1]] - lods[[progeny]][[2]])
  p <- ggplot(d) + geom_tile(aes(x=Var1, y=Var2, fill=value))
  p <- p + scale_fill_gradient2("log odds QQ/UU") + xlab("parent 1") + ylab("parent 2")
  p
})
