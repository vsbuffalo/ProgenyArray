# utilities.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

extractCutLabels <- function(x) {
    read.csv(text = gsub("\\[|\\]|\\(|\\)", "", x), header=FALSE,
                    col.names=c("start", "end"))
}



#' Create a genotyping error matrix.
#'
#' Returns a matrix of transition probabilities for a given genotype given error.
#'
#' @param ehet the probability that a heterozygous genotype is incorrect
#' @param ehom the probability that a homozygous genotype is incorrect
#'
#' This function creates a matrix of the form:
#'   [1-e, e/2, e/2,
#'    E/2, 1-E, E/2,
#'    e/2, e/2, 1-e]
#' where E is the probability the heterozygous genotype call is incorrect, and
#' e is the probability that the homozygous genotypes call is incorrect.
#' @export
genotypingErrorMatrix <- function(ehet=0.6, ehom=0.1) {
  matrix(c(1-ehom, ehet/2, ehom/2, ehom/2, 1-ehet, ehom/2, ehom/2, ehet/2, 1-ehom),
         ncol=3)
}

#' Calculate allele frequenceies.
#'
#' @param x a genotype matrix
#' @param min force minimum allele frequency to be 1/nsample
#' @export
alleleFreqs <- function(x, min=FALSE) {
  # TODO make this handle missing data better.
  nind <- ncol(x)
  nmissing <- rowSums(is.na(x))
  out <- rowSums(x, na.rm=TRUE)/(2*(nind - nmissing))
  if (min)
    return(ifelse(out == 0, 1/ncol(x), out))
  return(out)
}


#' Calculate allele frequencies from parents in a ProgenyArray
#'
#' @param x a ProgenyArray object
#'
#' @export
setMethod("calcFreqs", "ProgenyArray",
          function(x) {
            x@freqs <- alleleFreqs(parentGenotypes(x))
            x
})

#' Getter for allele frequencies
#'
#' @param x a ProgenyArray object
setMethod("freqs", "ProgenyArray",
          function(x) {
            if (length(freqs) == 0)
              stop("allele frequencies have not been calculated, run calcFreqs")
            x@freqs
})

#' Calculate genotype frequencies.
#'
#' @param x a genotype matrix
#' @param counts whether to return genotype counts rather than frequency
#' @export
genotypeFreqs <- function(x, counts=FALSE) {
  out <- t(apply(x, 1, countGenotypes))
  if (!counts)
    return(sweep(out, 1, rowSums(tt), "/"))
  return(out)
}

#' Create dataframe of Hardy-Weinberg values.
#'
#' @param x a genotype matrix
#' @export
HW <- function(x) {
  out <- as.data.frame(genotypeFreqs(x))
  colnames(out) <- c("g0", "g1", "g2", "gNA")
  out$freq <- alleleFreqs(x)
  out
}

#' Plot genotype and allele frequencies.
#'
#' @param x a genotype matrix
#' @export
HWPlot <- function(x) {
  p <- ggplot(HW(x))
  p <- p + geom_point(aes(x=freq, y=g0), color="red")
  p <- p + geom_point(aes(x=freq, y=g1), color="blue")
  p <- p + geom_point(aes(x=freq, y=g2), color="dark green")
  p
}

MendelianTransmissionList <- function(x) {
    # Mendelian transition matrix, in list notation for readability
    # TODO: we can build this up programmatically using Kronecker products.
    # first dimension (list element): mom's genotype
    # second dimension: alleged dad's genotype
    # third dimension: progeny genotype
    MT_list <- list(
               "0"=matrix(c(1, 1/2, 0, 0, 1/2, 1, 0, 0, 0),
                          ncol=3),
               "1"=matrix(c(1/2, 1/4, 0, 1/2, 1/2, 1/2, 0, 1/4, 1/2),
                          ncol=3),
               "2"=matrix(c(0, 0, 0, 1, 1/2, 0, 0, 1/2, 1),
                          ncol=3))
    return(MT_list)
}

MendelianTransmissionMatrix <- function() {
  # convert to array. Now dimensions are:
  # (1) father, (2) progeny, (3) mother
  m <- MendelianTransmissionList()
  MT <- array(unlist(m), dim=c(3, 3, 3))
  stopifnot(all(sapply(1:3, function(i) m[[i]] == MT[, , i])))
  return(MT)
}

validGenotypes <- function(x, allowNA=TRUE) {
  genos <- seq(0, 2)
  if (allowNA)
    genos <- c(NA, genos)
  all(x %in% genos)
}

MendelianInconsistencies <- function(parent1, child, parent2=NULL) {
  m <- MendelianTransmissionMatrix() > 0
  stopifnot(validGenotypes(c(parent1, child)))
  parent1 <- parent1+1L
  child <- child+1L
  if (is.null(parent2)) {
    validGenotypes(parent2)
    return(mapply(function(p1, c) all(!m[, c, p1]), parent1, child))
  }
  parent2 <- parent2+1L
  return(mapply(function(p1, p2, c) !m[p2, c, p1], parent1, parent2, child))
}

#' Trio Mendelian inconsistency across full-sib family per phased tile
#'
#' @param parent_1 parent 1
#' @param parent_2 parent 2
#' @param chrom which chromosome was used
#' @param tiles tiles used in phasing
#' @param phased_1 phasing results for parent 1 (if \code{NULL}, uses genotype data)
#' @param phased_2 phasing results for parent 2 (if \code{NULL}, uses genotype data)
#' @export
miMatrixPerTile <- function(x, parent_1, parent_2, chrom, tiles, phased_1=NULL, phased_2=NULL) {
  mi <- lapply(1:length(tiles), function(i) {
    tile_i <- tiles[[i]]
    if (!is.null(phased_1)) {
      p1_geno <- rowSums(phased_1[[i]]$haplos)
    } else {
      p1_geno <- unname(parentGenotypes(x, seqname=chrom)[tile_i, parent_1])
    }
    if (!is.null(phased_2)){
      p2_geno <- rowSums(phased_2[[i]]$haplos)
    } else {
      p2_geno <- unname(parentGenotypes(x, seqname=chrom)[tile_i, parent_2])
    }
    pars <- parents(x)
    progi <- pars$progeny[pars$parent_1 == parent_1 & pars$parent_2 == parent_2]
    pgeno <- progenyGenotypes(x, seqname=chrom)[tile_i, progi]
    stopifnot(nrow(pgeno) == length(p1_geno))
    stopifnot(length(p2_geno) == length(p1_geno))
    apply(pgeno, 2, function(child) MendelianInconsistencies(p1_geno, child, p2_geno))
  })
  mi
}

#' Mendelian inconsistency rates across a chromosome
#'
#' @param mi results from \code{miMatrixPerTile}
#' @param x ProgenyArray object
#' @param chrom which chromosome was used
#' @param tiles tiles used in phasing
miRates <- function(mi, x, chrom, tiles) {
  tmp <- lapply(1:length(mi), function(i) {
                m <- mi[[i]]
                tmp <- setNames(melt(m), c("index", "progeny", "value"))
                out <- summarise(group_by(tmp, 'index'),
                          inconsistent=sum(value, na.rm=TRUE)/sum(!is.na(value)),
                          na=sum(is.na(value))/length(value))
                out$tile <- i
                out
  })
  d <- do.call(rbind, tmp)
  d$pos <- start(x@ranges[as.logical(seqnames(x@ranges) == chrom)])
  d
}


#### BEGIN DEPRECATED ####
# These functions were used to generate Mendelian transmission probabilties
# when parentage was inferred in R. Now this is done in C++; these are just
# here for debugging/testing


conditionalOffspringParentMatrix <- function(freq) {
  #  P(Go | Gm)
  #      offspring
  #    0   1          2
  # 0  b   c          0
  # 1  b/2  (b+c)/2   c/2
  # 2  0    b         c
  b <- 1 - freq
  c <- freq
  matrix(c(b, b/2, 0, c, (b+c)/2, b, 0, c/2, c), nrow=3)
}


#' Create a Mendelian transmission matrix, with error in observed progeny genotypes.
#'
#' @param ehet the probability that a heterozygous genotype is incorrect
#' @param ehom the probability that a homozygous genotype is incorrect
#' @export
MendelianTransmissionWithError <- function(ehet=0.6, ehom=0.1) {
  out <- lapply(MendelianTransmissionList(),
                function(x) x %*% genotypingErrorMatrix(ehet, ehom))
  array(unlist(out), dim=c(3, 3, 3))
}

probOffspringGivenParent <- function(offspring, parent, freqs, ehet, ehom) {
  stopifnot(!any(is.na(offspring)))
  stopifnot(!any(is.na(parent)))
  stopifnot(!any(is.na(freqs)))
  error_matrix <- genotypingErrorMatrix(ehet, ehom)
  mapply(function(f, o, m) {
         mat <- conditionalOffspringParentMatrix(f)
         (mat %*% error_matrix)[m+1L, o+1L]
  }, freqs, offspring, parent)
}


#### END DEPRECATED ####
