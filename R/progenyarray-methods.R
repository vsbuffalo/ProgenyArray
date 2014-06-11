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
ProgenyArray <- function(geno, loci=GRanges(), ref=character(),
                         alt=CharacterList(), samples=character()) {
	obj <- new("ProgenyArray", genotypes=geno, ranges=loci,
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
			ref=tasselr:::ref(from), alt=tasselr:::alt(from), samples=from@samples)
})


#' Pretty-print a ProgenyArray object
#'
#' @param object a ProgenyArray object
#'
#' @export
setMethod("show",
          c(object="ProgenyArray"),
          function(object) {
            cat(sprintf("ProgenyArray object\n%d loci x %d samples\n",
                        nrow(object@genotypes), ncol(object@genotypes)))
            cat(sprintf("Number of chromosomes: %d\nObject size: %s Mb\n",
                        length(seqlevels(object@ranges)),
                        round(object.size(object)/1024^2, 3)))
            numOrNA <- function(x) ifelse(length(x) == 0, NA, length(x))
            nfathers <- numOrNA(object@fathers)
            nposs_fathers <- numOrNA(object@possible_fathers)
            nmothers <- numOrNA(object@mothers)
            nprogeny <- numOrNA(object@progeny)
            cat(sprintf("Number of possible fathers: %d\n", nposs_fathers))
            cat(sprintf("Number of fathers: %d\n", nfathers))
            cat(sprintf("Number of mothers: %d\n", nmothers))
            cat(sprintf("Number of progeny: %d\n", nprogeny))
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

#' Accessor for progeny in a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @export
setMethod("progeny",
          c(x="ProgenyArray"),
          function(x) {
            return(x@progeny)
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

#' Set method for fathers
#'
#' @param object
#' @param value
#' @export
#' @name set-methods
setReplaceMethod("fathers", "ProgenyArray", function(object, value) {
	object@fathers <- value
	return(object)
})

#' Set method for possibleFathers
#'
#' @param object
#' @param value
#' @export
#' @name set-methods
setReplaceMethod("possibleFathers", "ProgenyArray", function(object, value) {
	object@possible_fathers <- value
	return(object)
})


#' Set method for mothers
#'
#' @param object
#' @param value
#' @export
#' @name set-methods
setReplaceMethod("mothers", "ProgenyArray", function(object, value) {
	object@mothers <- value
	return(object)
})


#' Set method for progeny
#'
#' @param object
#' @param value
#' @export
#' @name set-methods
setReplaceMethod("progeny", "ProgenyArray", function(object, value) {
	object@progeny <- value
	return(object)
})

