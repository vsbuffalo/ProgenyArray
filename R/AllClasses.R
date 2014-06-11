# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' An S4 class that stores genotyping data from a progeny array design.
#'
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref reference alleles
#' @slot alt alternate alleles
#' @slot genotypes a matrix of bialleic genotypes
#' @slot progeny an integer vector indicating which samples are progeny
#' @slot possible_fathers an integer vector indicating which samples are possible fathers
#' @slot mothers integer vector indicating the mother of each progeny
#' @slot fathers integer vector indicating the father of each progeny
#' @slot samples sample names
#'
#' @exportClass ProgenyArray
setClass("ProgenyArray",
         representation=representation(ranges="GRanges",
                                       ref="character",
                                       alt="CharacterList",
                                       genotypes="matrix",
                                       progeny="integer",
                                       possible_fathers="integer",
                                       mothers="integer",
                                       fathers="integer",
                                       fathers_lle="list",
                                       samples="character"),
         prototype=prototype(ranges=GRanges(),
                             ref=character(),
                             alt=CharacterList(),
                             progeny=integer(),
                             possible_fathers=integer(),
                             mothers=integer(),
                             fathers=integer(),
                             fathers_lle=list(),
                             samples=character()),
         validity=function(object) {
           if (length(mothers(object) != length(progeny)))
             stop("mothers vector must be same length as progeny vector")
         })

