# AllGenerics.R --

# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' Method to calculate if progeny is selfed
#'
setGeneric("isSelfed", function(x) {
           standardGeneric("isSelfed")
             })

#' Accessor for alternate alleles
#'
setGeneric("alt", function(x) {
           standardGeneric("alt")
             })

#' Accessor for reference alleles
#'
setGeneric("ref", function(x) {
           standardGeneric("ref")
             })

#' Accessor for fathers (returning indices or names)
#'
setGeneric("fathers", function(x, names=TRUE) {
           standardGeneric("fathers")
             })

#' Accessor for mothers (returning indices or names)
#'
setGeneric("mothers", function(x, names=TRUE) {
           standardGeneric("mothers")
             })

#' Accessor for user-supplied mothers (returning indices or names)
#'
setGeneric("suppliedMothers", function(x) {
           standardGeneric("suppliedMothers")
             })

#' Accessor for parents (returning indices or names)
#'
setGeneric("parents", function(x, names=TRUE) {
           standardGeneric("parents")
             })

#' Accessor for progenyGenotypes
#'
setGeneric("progenyGenotypes", function(x, seqname=NULL) {
           standardGeneric("progenyGenotypes")
             })

#' Accessor for parentGenotypes
#'
setGeneric("parentGenotypes", function(x, seqname=NULL) {
           standardGeneric("parentGenotypes")
             })

#' Accessor for inferParents
#'
setGeneric("inferParents", function(x, ehet, ehom, verbose=FALSE) {
           standardGeneric("inferParents")
             })

#' Accessor for progeny sample names
#'
setGeneric("progenyNames", function(x) {
           standardGeneric("progenyNames")
             })

#' Accessor for parents sample names
#'
setGeneric("parentNames", function(x) {
           standardGeneric("parentNames")
             })

#' Create a diagnostic plot of a parentage inference, looking at all log likelihood ratios.
#'
setGeneric("parentLLPlot", function(x, progeny) {
           standardGeneric("parentLLPlot")
             })

#' Setter for progeny sample names
#'
setGeneric("progenyNames<-", function(x, value) {
           standardGeneric("progenyNames<-")
             })

#' Setter for parents sample names
#'
setGeneric("parentNames<-", function(x, value) {
           standardGeneric("parentNames<-")
             })

#' Accessor for user-supplied mothers
#'
setGeneric("suppliedMothers<-", function(x, value) {
           standardGeneric("suppliedMothers<-")
             })

#' Setter for mothers
#'
setGeneric("mothers<-", function(x, value) {
           standardGeneric("mothers<-")
             })

#' Setter for fathers
#'
setGeneric("fathers<-", function(x, value) {
           standardGeneric("fathers<-")
             })

#' Create parentage dataframe
#'
setGeneric("parentage", function(x) {
           standardGeneric("parentage")
             })

#' Setter for LOD cutoff
#'
setGeneric("lodcutoff<-", function(x, value) {
           standardGeneric("lodcutoff<-")
             })

#' Getter for LOD cutoff
#'
setGeneric("lodcutoff", function(x) {
           standardGeneric("lodcutoff")
             })

#' Calculate allele frequencies
#'
setGeneric("calcFreqs", function(x) {
           standardGeneric("calcFreqs")
             })

#' Getter for frequencies
#'
setGeneric("freqs", function(x) {
           standardGeneric("freqs")
             })

#' Setter for frequencies
#'
setGeneric("freqs<-", function(x, value) {
           standardGeneric("freqs<-")
             })

#' Phase all parents data in a ProgenyArray
#'
setGeneric("phaseParents", function(x, tiles, ehet=0.8, ehom=0.1,
                                    na_thresh=0.8, parallel=FALSE, verbose=TRUE) {
    standardGeneric("phaseParents")
})

#### Generics for SimulationData

#' Getter for ProgenyArray object from a SimulationData object
setGeneric("progenyArray", function(x) {
           standardGeneric("progenyArray")
          })


