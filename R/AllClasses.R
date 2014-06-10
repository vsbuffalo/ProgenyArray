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
         slots=list(ranges="GRanges",
                    ref="character",
                    alt="CharacterList",
                    genotypes="matrix",
										progeny="integer",
										possible_fathers="integer",
										mothers="integer",
										fathers="integer",
                    samples="character"))
