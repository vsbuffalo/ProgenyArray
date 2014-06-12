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
#' @slot mothers integer vector indicating the mother of each progeny
#' @slot fathers integer vector indicating the father of each progeny
#' @slot progeny_samples sample names of progeny
#' @slot parents_samples sample names of parents
#' @slot complete_loci which loci are complete in the parents
#'
#' @exportClass ProgenyArray
setClass("ProgenyArray",
         representation=representation(ranges="GRanges",
                                       ref="character",
                                       alt="CharacterList",
                                       progeny_geno="matrix",
                                       parents_geno="matrix",
                                       mothers="integer",
                                       fathers="integer",
                                       fathers_lle="list",
                                       progeny_samples="character",
                                       parents_samples="character",
																			 complete_loci="integer"),
         prototype=prototype(ranges=GRanges(),
                             ref=character(),
                             alt=CharacterList(),
                             mothers=integer(),
                             fathers=integer(),
                             fathers_lle=list(),
                             progeny_samples=character(),
                             parents_samples=character(),
														 complete_loci=integer()),
         validity=function(object) {
           if (nrow(object@progeny_geno) != nrow(object@parents_geno))
             stop("number of progeny loci must equal number of parent loci")
           if (length(object@mothers) != ncol(object@progeny_geno))
             stop("mothers vector must be equal length as progeny genotype matrix")
         })

