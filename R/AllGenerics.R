# AllGenerics.R -- 
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.


#' Accessor for alternate alleles
#'
setGeneric("alt", function(x, as_char=TRUE) {
           standardGeneric("alt")
             })

#' Accessor for reference alleles
#'
setGeneric("ref", function(x, as_char=TRUE) {
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
setGeneric("suppliedMothers", function(x, names=TRUE) {
           standardGeneric("suppliedMothers")
             })

#' Accessor for parents (returning indices or names)
#'
setGeneric("parents", function(x, names=TRUE) {
           standardGeneric("parents")
             })


#' Accessor for progenyGenotypes
#'
setGeneric("progenyGenotypes", function(x, names=TRUE) {
           standardGeneric("progenyGenotypes")
             })

#' Accessor for parentGenotypes
#'
setGeneric("parentGenotypes", function(x, names=TRUE) {
           standardGeneric("parentGenotypes")
             })

#' Accessor for inferParents
#'
setGeneric("inferParents", function(x, ehet, ehom, verbose=FALSE) {
           standardGeneric("inferParents")
             })

#' Accessor for progeny sample names
#'
setGeneric("progenyNames", function(object, value) {
           standardGeneric("progenyNames")
             })

#' Accessor for parents sample names
#'
setGeneric("parentNames", function(object, value) {
           standardGeneric("parentNames")
             })

#' Create a diagnostic plot of a parentage inference, looking at all log likelihood ratios.
#'
setGeneric("parentLLPlot", function(x, progeny) {
           standardGeneric("parentLLPlot")
             })

#' Setter for progeny sample names
#'
setGeneric("progenyNames<-", function(object, value) {
           standardGeneric("progenyNames<-")
             })

#' Setter for parents sample names
#'
setGeneric("parentNames<-", function(object, value) {
           standardGeneric("parentNames<-")
             })

#' Accessor for user-supplied mothers
#'
setGeneric("suppliedMothers<-", function(object, value) {
           standardGeneric("suppliedMothers<-")
             })

#' Accessor for mothers
#'
setGeneric("mothers<-", function(object, value) {
           standardGeneric("mothers<-")
             })

#' Accessor for fathers
#'
setGeneric("fathers<-", function(object, value) {
           standardGeneric("fathers<-")
             })

#' Create parentage dataframe
#'
setGeneric("parentage", function(x) {
           standardGeneric("parentage")
             })

