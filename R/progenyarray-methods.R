# progenyarray-methods.R -- method for working with ProgenyArray objects
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' Return which loci are complete
whichLociComplete <- function(x) {
  which(apply(!is.na(x), 1, all))
}

#' Constructor new ProgenyArray object.
#'
#' Create a new ProgenyArray object for storing genotyping and other data from
#' a progeny array experiment.
#'
#' @param progeny_geno
#' @param parents_geno
#' @param loci a \code{GRanges} object of loci
#' @param ref a character vector of reference loci
#' @param alt a \code{CharacterList} of alternate loci
#' @export
ProgenyArray <- function(progeny_geno, parents_geno, mothers=integer(),
                         loci=GRanges(), ref=character(), alt=CharacterList()) {
  obj <- new("ProgenyArray", progeny_geno=progeny_geno, parents_geno=parents_geno,
             mothers=mothers, ranges=loci, ref=ref, alt=alt)

  obj@complete_loci <- whichLociComplete(obj@parents_geno)
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
            cat(sprintf("Number of chromosomes: %d\nObject size: %s Mb\n",
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
            cat(sprintf("Number of progeny: %d\n", nprogeny))
            cat(sprintf("Number of parents: %d\n", nparents))
            #cat(sprintf("Number of fathers: %d\n", nfathers))
            #cat(sprintf("Number of mothers: %d\n", nmothers))
            cat(sprintf("Proportion missing:\n  progeny: %0.3f\n  parents: %0.3f\n",
                        sum(is.na(object@progeny_geno))/length(object@progeny_geno),
                        sum(is.na(object@parents_geno))/length(object@parents_geno)))
            cat(sprintf("Number of complete parental loci: %d\n",
                        length(object@complete_loci)))

          })


#' Accessor for parent genotypes
#'
#' @param x a ProgenyArray object
#' @export
setMethod("parentGenotypes",
          c(x="ProgenyArray"),
          function(x) {
            return(x@parents_geno)
          })

#' Accessor for progeny genotypes
#'
#' @param x a ProgenyArray object
#' @export
setMethod("progenyGenotypes",
          c(x="ProgenyArray"),
          function(x) {
            return(x@progeny_geno)
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
#' @param as_char logical indicating whether to return as \code{CharacterList}
#' or as Tassel encodes them, as an \code{IntegerList}
#' @export
setMethod("alt",
          c(x="ProgenyArray"),
          function(x, as_char=TRUE) {
            if (!as_char)
              return(x@alt)
            out <- lapply(x@alt, function(x) TASSELL_ALLELES[x])
            return(as(out, "CharacterList"))
          })

#' Accessor for reference alleles from a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @param as_char logical indicating whether to return as character
#' or as Tassel encodes them, as integer
#' @export
setMethod("ref",
          c(x="ProgenyArray"),
          function(x, as_char=TRUE) {
            if (!as_char)
              return(x@ref)
            return(TASSELL_ALLELES[x@ref])
          })

#' Return progeny sample names
#'
#' @export
setMethod("progenyNames", "ProgenyArray", function(object) {
  colnames(object@progeny_geno)
})

#' Return parent sample names
#'
#' @export
setMethod("parentNames", "ProgenyArray", function(object) {
  colnames(object@parents_geno)
})

##
###' Set method for fathers
###'
###' @name fathers
###' @export
##setReplaceMethod("fathers", "ProgenyArray", function(object, value) {
##  if (length(value) != ncol(object@progeny))
##    stop("length of value must be same as number of progeny")
##  object@fathers <- value
##  return(object)
##})
##
###' Set method for mothers
###'
###' @name mothers
###' @export
##setReplaceMethod("mothers", "ProgenyArray", function(object, value) {
##  if (length(value) != ncol(object@progeny_geno))
##    stop("length of value must be same as number of progeny")
##  object@mothers <- value
##  return(object)
##})

#' Set method for user suppied-mothers
#'
#' @name mothers
#' @export
setReplaceMethod("suppliedMothers", "ProgenyArray", function(object, value) {
  if (length(value) != ncol(object@progeny_geno))
    stop("length of value must be same as number of progeny")
  object@supplied_mothers <- value
  return(object)
})


#' Set method for parents sample names
#'
#' @name parentNames
#' @export
setReplaceMethod("parentNames", "ProgenyArray", function(object, value) {
  if (length(value) != ncol(object@parents_geno))
    stop("length of value must be same as number of parents")
  if (any(duplicated(value)))
    stop("no duplicated names allowed")
  colnames(object@parents_geno) <- value
  return(object)
})


#' Set method for progeny sample names
#'
#' @name progenyNames
#' @export
setReplaceMethod("progenyNames", "ProgenyArray", function(object, value) {
  if (length(value) != ncol(object@progeny_geno))
    stop("length of value must be same as number of progeny")
  if (any(duplicated(value)))
    stop("no duplicated names allowed")
  colnames(object@progeny_geno) <- value
  return(object)
})

