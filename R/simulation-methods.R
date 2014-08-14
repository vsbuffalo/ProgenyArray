# simulation-methods.R -- methods for working with the SimulationData objects
# these simulation methods do not keep track of haplotypes, and are just
# appropriate for parentage methods.

# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

sampleSiteFreqs <- function(n, nsites) {
	sf <- 1/seq_len(n) # sfs for sample of n
	site_probs <- sf/sum(sf) # probability of getting a certain number of sites
	site_freqs <- seq_len(n)/n
	sample(site_freqs, nsites, replace=TRUE, prob=site_probs)
}

createDiploidSample  <- function(n, nsites) {
	# this is a bit of an approximation: create SFS, and then generate site
	# frequencies from that.
  sample_sf <- sampleSiteFreqs(n, nsites)

  #samples <- t(replicate(2*n, { Vectorize(function(prob) rbinom(1, 1, prob), "prob")(sample_sf) } ))
	#split(split(samples, row(samples)), rep(seq_len(n), each=2))
  out <- lapply(seq_len(n), function(x) rbinom(nsites, 1, sample_sf)+rbinom(nsites, 1, sample_sf))
  do.call(cbind, out)
}

mate <- function(x, y, dont_combine=FALSE) {
  # mate two genotype vectors, each genotype is independently inherited
  transmit <- list(c(0), c(0, 1), c(1, 1))
  x_gametes <- sapply(transmit[x+1L], sample, size=1)
  y_gametes <- sapply(transmit[y+1L], sample, size=1)
  if (!dont_combine)
    return(x_gametes + y_gametes)
  cbind(x_gametes, y_gametes)
}

addGenotypeError <- function (x, ehet=0, ehom=0) {
  genotype_error <- genotypingErrorMatrix(ehet, ehom)
  sapply(x, function(y) sample(c(0L, 1L, 2L), 1, prob=genotype_error[y+1L, ]))
}

createProgenyArrayWithSelfing <- function(x, nprogeny=50, selfing=0.5) {
  # create progeny array for a single mother from a population x
  inds <- seq_len(ncol(x))
  mom_i <- sample(inds, nprogeny, replace=TRUE)
  isselfed <- as.logical(rbinom(nprogeny, 1, selfing))
  dad_i <- mapply(function(s, m) {
                  if (s) # dad is mom
                    return(m)
                  dad <- sample(setdiff(inds, m), 1, replace=TRUE) # outcrossed, exclude mom
                  return(dad)
                 }, isselfed, mom_i)
  # no accidental selfs
  stopifnot(all((mom_i == dad_i) == isselfed))
  # mate both parents using the mate function
  progeny <- do.call(cbind, mapply(function(m, d) mate(x[, m], x[, d]), mom_i, dad_i, SIMPLIFY=FALSE))
  list(progeny=progeny, mom_i=mom_i, dad_i=dad_i)
}

addMissingness <- function(x, rate) {
  ifelse(rbinom(length(x), 1, rate), NA, x)
}

#' Create some simulation data
#'
#' @param nparent number of parents
#' @param nprogeny number of progeny
#' @param nloci number of loci
#' @param selfing selfing rate
#' @param prop_parent_missing proportion of random missing data in parents
#' @param prop_progeny_missing proportion of random missing data in progeny
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#'
#' @export
SimulatedProgenyArray <- function(nparent, nprogeny, nloci, selfing=0.5,
                                  prop_parent_missing=0.2, prop_progeny_missing=0.2,
                                  ehet=0.8, ehom=0.1) {
  pop <- createDiploidSample(nparent, nloci)
  prog <- createProgenyArrayWithSelfing(pop, nprogeny, selfing)
  progeny_geno <- apply(prog$progeny, 2, function(x) addGenotypeError(x, ehet, ehom))
  parent_geno <- apply(pop, 2, function(x) addGenotypeError(x, ehet, ehom))
  progeny_genoNA <- apply(progeny_geno, 2, function(x) addMissingness(x, prop_progeny_missing))
  parent_genoNA <- apply(parent_geno, 2, function(x) addMissingness(x, prop_parent_missing))
  pa <- ProgenyArray(progeny_genoNA, parent_genoNA, mothers=prog$mom_i)
  new("SimulatedProgenyArray", nparent=as.integer(nparent), nprogeny=as.integer(nprogeny),
      nloci=as.integer(nloci),
      selfing=selfing,
      prop_parent_missing=prop_parent_missing,
      prop_progeny_missing=prop_progeny_missing,
      ehet=ehet, ehom=ehom,
      parent_geno=pop,
      progeny_geno=prog$progeny,
      mothers=prog$mom_i,
      fathers=prog$dad_i,
      progeny_array=pa)
}

#' Infer parents for SimulatedProgenyArray
#'
#' @param x a SimulatedProgenyArray object
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
setMethod("inferParents", "SimulatedProgenyArray", function(x, ehet, ehom, verbose=FALSE) {
  x@progeny_array <- inferParents(x@progeny_array, ehet, ehom, verbose)
  x
})

numCorrectParents <- function(real, inferred) {
  mapply(function(x, y) length(intersect(x, y)), real, inferred)
}

correctParentCounts <- function(real, inferred) {
  ff <- factor(numCorrectParents(real, inferred),
               levels=c(0, 1, 2))
  tbl <- table(ff)
  tbl
}

propCorrectParents <- function(x) {
  n <- ncol(x@progeny_geno)
  real <- split(cbind(x@mothers, x@fathers), seq_len(n))
  inferred <- split(x@progeny_array@parents[, 2:3], seq_len(n))
  tbl <- prop.table(correctParentCounts(real, inferred))
  unname(tbl[3])
}

setMethod("show", "SimulatedProgenyArray", function(object) {
  o <- object
  cat(sprintf("SimulatedProgenyArray object: %d loci, %d parents, %d progeny\n",
              o@nloci, o@nparent, o@nprogeny))
  cat("Targetted rates:\n")
  cat(sprintf("  selfing: %0.3f\n", o@selfing))
  cat(sprintf("  parent missing: %0.3f\n", o@prop_parent_missing))
  cat(sprintf("  progeny missing: %0.3f\n", o@prop_progeny_missing))
  cat(sprintf("  heterozygous error: %0.3f\n", o@ehet))
  cat(sprintf("  homozygous error: %0.3f\n", o@ehom))
  cat("*** Internal ProgenyArray Object ***\n")
  show(o@progeny_array)
  cat("*** Parentage Accuracy ***\n")
  if (nrow(progenyArray(o)@parents) > 0) {
    #mc <- sum(o@mothers == mothers(progenyArray(o)))
    #fc <- sum(o@fathers == fathers(progenyArray(o)))
    n <- ncol(o@progeny_geno)
    #cat(sprintf("%d/%d (%0.0f%%) mothers correct\n", mc, n, 100*mc/n))
    #cat(sprintf("%d/%d (%0.0f%%) fathers correct\n", fc, n, 100*fc/n))
    real <- split(cbind(o@mothers, o@fathers), seq_len(n))
    inferred <- split(o@progeny_array@parents[, 2:3], seq_len(n))
    tbl <- prop.table(correctParentCounts(real, inferred))
    cat(sprintf("prop parents correct:\n"))
    cat(sprintf("  0    1    2\n"))
    cat(sprintf("%0.2f %0.2f %0.02f\n", tbl[1], tbl[2], tbl[3]))
  }
})


#' Accessor for parent genotypes
#'
#' @param x a SimulationData object
#' @export
setMethod("parentGenotypes",
          c(x="SimulatedProgenyArray"),
          function(x) {
            return(x@parents_geno)
          })

#' Accessor for progeny genotypes
#'
#' @param x a SimulationData object
#' @export
setMethod("progenyGenotypes",
          c(x="SimulatedProgenyArray"),
          function(x) {
            return(x@progeny_geno)
          })

#' Accessor for ProgenyArray object
#'
#' @param x a SimulationData object
#' @export
setMethod("progenyArray",
          c(x="SimulatedProgenyArray"),
          function(x) {
            return(x@progeny_array)
          })

