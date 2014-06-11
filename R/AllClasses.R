# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' An S4 class that stores genotyping data from a progeny array design.
#'
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref reference alleles
#' @slot alt alternate alleles
#' @slot progeny a matrix of bialleic progeny genotypes
#' @slot parents a matrix of bialleic parents genotypes
#' @slot mothers integer vector indicating the mother of each progeny
#' @slot fathers integer vector indicating the father of each progeny
#' @slot samples sample names
#'
#' @exportClass ProgenyArray
setClass("ProgenyArray",
         representation=representation(ranges="GRanges",
                                       ref="character",
                                       alt="CharacterList",
                                       progeny="matrix",
                                       parents="matrix",
                                       mothers="integer",
                                       fathers="integer",
                                       fathers_lle="list",
                                       samples="character"),
         prototype=prototype(ranges=GRanges(),
                             ref=character(),
                             alt=CharacterList(),
                             mothers=integer(),
                             fathers=integer(),
                             fathers_lle=list(),
                             samples=character()),
         validity=function(object) {
           if (nrow(object@progeny) != nrow(object@parents))
             stop("inconsistent number of loci: nrow(progeny) != nrow(parents)")
           if (length(object@mothers) != ncol(object@progeny))
             stop("mothers vector must be same length as progeny vector")
         })

