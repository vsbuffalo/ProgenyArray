# progenyarray-methods.R -- method for working with ProgenyArray objects
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' Constructor new ProgenyArray object.
#'
#' Create a new ProgenyArray object for storing genotyping and other data from
#' a progeny array experiment.
#'
#' @param geno a genotype matrix
#' @param loci a \S4Class{GRanges} object of loci
#' @param ref a character vector of reference loci
#' @param alt a \S4Class{CharacterList} of alternate loci
#' @param a character vector of sample names
ProgenyArray <- function(geno, loci, ref, alt, samples) {
	obj <- new("ProgenyArray", genotypes=gneo, ranges=loci,
						 ref=ref, alt=alt, samples=samples)
	obj
}

#' Coercion method to create a ProgenyArray object from GBS data from Tassel.
#' @name as
#'
#' @param from a TasselHDF5 object
#' @export
#'
setAs("TasselHDF5", "ProgenyArray", function(from) {
	new("ProgenyArray", genotypes=from@genotypes, ranges=from@ranges,
			ref=ref(x), alt=alt(x), samples=samples)
})

#' Accessor for possible fathers in a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("possibleFathers",
          c(x="ProgenyArray"),
          function(x) {
            return(x@possible_fathers)
          })


#' Accessor for fathers in a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("fathers",
          c(x="ProgenyArray"),
          function(x) {
            return(x@fathers)
          })

#' Accessor for mothers in a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("mothers",
          c(x="ProgenyArray"),
          function(x) {
            return(x@mothers)
          })

#' Accessor for genotype information from a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("geno",
          c(x="ProgenyArray"),
          function(x) {
            return(x@genotypes)
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


#' Accessor for samples from a ProgenyArray object
#'
#' @param object a ProgenyArray object
#' @export
setMethod("samples",
          c(object="ProgenyArray"),
          function(object) {
						return(object@samples)
          })


