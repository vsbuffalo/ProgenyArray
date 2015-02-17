# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' A container for tiles used for phasing by imputation
#'
#' @slot type a character describing the type of tile (position, SNP, genetic)
#' @slot width the width of each tile
#' @slot tiles a list of lists, with the first list containing all chromosomes' tiles, and the second list containing the indices of the SNPs per tile
#'  @slot info a list containing supplementary info about tiles, e.g genetic map and related smoothed versions
#' @exportClass PhasingTiles
setClass("PhasingTiles",
         representation=representation(type="character",
             width="integer",
             tiles="list",
             info="list"),             
             prototype=prototype(type=character(),
                 width=integer(),
                 tiles=list(), info=list()))

#' An S4 class that stores genotyping data from a progeny array design.
#'
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref reference alleles
#' @slot alt alternate alleles
#' @slot progeny_geno a matrix of bialleic progeny genotypes
#' @slot parents_geno a matrix of bialleic parents genotypes
#' @slot supplied_mothers integer vector indicating the mother of each progeny, supplied by user
#' @slot mothers integer vector indicating the inferred mother of each progeny
#' @slot fathers integer vector indicating the inferred father of each progeny
#' @slot parents dataframe with inferred parents and model information
#' @slot complete_loci which loci are not all missing or fixed in the parents
#' @slot parent_lods matrices from C++ parentage routine (for debugging)
#' @slot lodcutoff the LOD cutoff used in parentage inference
#' @slot phased_parents a list of parents' phased tiles for all chromosomes
#' @slot phased_parents_metadata a list of phasing parameters
#' @slot tiles a list of tiles for all chromosomes
#'
#' @exportClass ProgenyArray
setClass("ProgenyArray",
         representation=representation(ranges="GRanges",
             ref="character",
             alt="character",
             progeny_geno="matrix",
             parents_geno="matrix",
             supplied_mothers="integer", # TODO we can get rid of this
             mothers="integer",
             fathers="integer",
             freqs="numeric",
             parents="data.frame",
             complete_loci="integer",
             parent_lods="list",
             lodcutoff="numeric",
             # list of lists of lists (parents, chroms, and tiles)
             phased_parents="list",
             phased_parents_metadata="list",
             tiles="PhasingTiles",
             ligation="integer"),
         prototype=prototype(ranges=GRanges(),
             ref=character(),
             alt=character(),
             supplied_mothers=integer(),
             mothers=integer(),
             fathers=integer(),
             freqs=numeric(),
             parents=data.frame(),
             lodcutoff=NA_real_,
             phased_parents=list(),
             phased_parents_metadata=list(),
             tiles=NULL,
             ligation=NA_integer_),
         validity=function(object) {
             if (nrow(object@progeny_geno) != nrow(object@parents_geno))
                 stop("number of progeny loci must equal number of parent loci")
         })

#' An S4 class that wraps a ProgenyArray object of noisy, missing data
#'
#' @slot nparent parent population size
#' @slot nprogeny progeny population size
#' @slot nloci number of loci
#' @slot selfing selfing rate
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
