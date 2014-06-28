# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' An S4 class that stores genotyping data from a progeny array design.
#'
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref reference alleles
#' @slot alt alternate alleles
#' @slot progeny_geno a matrix of bialleic progeny genotypes
#' @slot parents_geno a matrix of bialleic parents genotypes
#' @slot supplied_mothers integer vector indicating the mother of each progeny, supplied by user
#' @slot mothers integer vector indicating the mother of each progeny
#' @slot fathers integer vector indicating the father of each progeny
#' @slot complete_loci which loci are complete in the parents
#'
#' @exportClass ProgenyArray
setClass("ProgenyArray",
         representation=representation(ranges="GRanges",
                                       ref="character",
                                       alt="CharacterList",
                                       progeny_geno="matrix",
                                       parents_geno="matrix",
                                       supplied_mothers="integer",
                                       mothers="integer",
                                       fathers="integer",
                                       parents="data.frame",
                                       fathers_lle="list",
                                       complete_loci="integer",
                                       parent_lods="list",
                                       lodcutoff="numeric"),
         prototype=prototype(ranges=GRanges(),
                             ref=character(),
                             alt=CharacterList(),
                             supplied_mothers=integer(),
                             mothers=integer(),
                             fathers=integer(),
                             parents=data.frame(),
                             fathers_lle=list(),
                             lodcutoff=NA_real_),
         validity=function(object) {
           if (nrow(object@progeny_geno) != nrow(object@parents_geno))
             stop("number of progeny loci must equal number of parent loci")
         })

