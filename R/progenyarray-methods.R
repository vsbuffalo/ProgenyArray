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
ProgenyArray <- function(progeny, parents, mothers=integer(), loci=GRanges(),
                         ref=character(), alt=CharacterList(),
                         samples=character()) {
	obj <- new("ProgenyArray", progeny=progeny, parents=parents,
             mothers=mothers, ranges=loci, ref=ref, alt=alt,
             samples=samples)
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
                        nrow(object@progeny), ncol(object@parents),
                        ncol(object@progeny)))
            cat(sprintf("Number of chromosomes: %d\nObject size: %s Mb\n",
                        length(seqlevels(object@ranges)),
                        round(object.size(object)/1024^2, 3)))
            numOrNA <- function(x) ifelse(length(x) == 0, NA, length(x))
            nfathers <- numOrNA(object@fathers)
            nparents <- numOrNA(object@parents)
            nmothers <- numOrNA(object@mothers)
            nprogeny <- numOrNA(object@progeny)
            cat(sprintf("Number of progeny: %d\n", nprogeny))
            cat(sprintf("Number of parents: %d\n", nparents))
            cat(sprintf("Number of fathers: %d\n", nfathers))
            cat(sprintf("Number of mothers: %d\n", nmothers))
          })

#' Accessor for parent genotypes
#'
#' @param x a ProgenyArray object
#' @export
setMethod("parentGenotypes",
          c(x="ProgenyArray"),
          function(x) {
            return(x@parents)
          })

#' Accessor for progeny genotypes
#'
#' @param x a ProgenyArray object
#' @export
setMethod("progenyGenotypes",
          c(x="ProgenyArray"),
          function(x) {
            return(x@progeny)
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

