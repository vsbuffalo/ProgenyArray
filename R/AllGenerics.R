# AllGenerics.R -- 
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.


#' Accessor for genotype matrix
#'
setGeneric("geno", function(x) {
           standardGeneric("geno")
             })

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

#' Accessor for progeny (returning indices or names)
#'
setGeneric("progeny", function(x, names=TRUE) {
           standardGeneric("progeny")
             })

#' Accessor for possibleFathers (returning indices or names)
#'
setGeneric("possibleFathers", function(x, names=TRUE) {
           standardGeneric("possibleFathers")
             })

#' Accessor for inferFathers
#'
setGeneric("inferFathers", function(x, ehet, ehom) {
           standardGeneric("inferFathers")
             })

#' Setter for progeny
#'
setGeneric("progeny<-", function(object, value) {
           standardGeneric("progeny<-")
             })

#' Setter for possibleFathers
#'
setGeneric("possibleFathers<-", function(object, value) {
           standardGeneric("possibleFathers<-")
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

