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

#' An S4 class that wraps a ProgenyArray object of noisy, missing data
#'
#' @slot nparent parent population size
#' @slot nprogeny progeny population size
#' @slot nloci number of loci
#' @slot prop_parent_missing proportion of genotype data that are missing in parents
#' @slot prop_progeny_missing proportion of genotype data that are missing in progeny
#' @slot ehet heterozygous error rate
#' @slot ehom homozygous error rate
#' @slot progeny_geno full genotype matrix for progeny
#' @slot parent_geno full genotype matrix for parents
#' @slot mothers true mothers of each progeny
#' @slot fathers true mothers of each progeny
#' @slot progeny_array ProgenyArray object of these data, with missing values and errors
#'
#' @exportClass SimulatedProgenyArray
setClass("SimulatedProgenyArray",
         representation=representation(nparent="integer",
                                       nprogeny="integer",
                                       nloci="integer",
                                       selfing="numeric",
                                       prop_parent_missing="numeric",
                                       prop_progeny_missing="numeric",
                                       ehet="numeric",
                                       ehom="numeric",
                                       # the full, non-missing genotype matrices
                                       parent_geno="matrix",
                                       progeny_geno="matrix",
                                       mothers="integer",
                                       fathers="integer",
                                       progeny_array="ProgenyArray"),
         prototype=prototype(nparent=integer(),
                             nprogeny=integer(),
                             nloci=integer(),
                             selfing=numeric(),
                             prop_parent_missing=numeric(),
                             prop_progeny_missing=numeric(),
                             ehet=numeric(),
                             ehom=numeric(),
                             # the full, non-missing genotype matrices
                             parent_geno=matrix(),
                             progeny_geno=matrix(),
                             mothers=integer(),
                             fathers=integer()))

