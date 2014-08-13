

sampleAllele <- function(x) {
  sapply(x, function(g) {
         if (g == 1) return(rbinom(1, 1, 0.5))
         if (g == 0) return(0L)
         if (g == 2) return(1L)
         stop("unknown genotype value specified; genotypes need to be in (0, 1, 2)")
  })
}

sampleHaplotypes <- function(m, f) {
  which <- setNames(sample(1:2, 2, replace=TRUE), c('m', 'f'))
  list(gametes=cbind(m[, which[1]], f[, which[2]]), which=which)
}

haplotypeSample2List <- function(x, label, ehet, ehom) {
  # convenience function that takes a list of individuals' gametes and converts
  # these to a genotype matrix (adding error) and a dataframe of which maternal
  # and paternal gametes were inherited.
  gametes <- lapply(x, `[[`, 1)
  d <- data.frame(do.call(rbind, lapply(x, `[[`, 2)))
  d$i <- seq_along(gametes)
  d$type <- label
  genos <- do.call(cbind, lapply(gametes,
             function(x) addGenotypeError(rowSums(x), ehet=ehet, ehom=ehom)))
  list(geno=genos, gametes=gametes, which=d)
}

#' Simulate a half-sib family, fixed number of sites
#'
#' @export n size of halfsib family
#' @export nsites number of sites
#' @export popsize population size (from which parents are sampled)
#' @export selfing selfing rate
#' @export fullsib fullsib rate (1 - (selfing + fullsib)) = halfsib rate
#' @export nfullsib number of fullsib families
#' @export mom_inbred whether mom is inbred
#' @export ehet heterozygous error rate
#' @export ehom homozygous error rate
#'
#' Returns a list with elements \code{geno}, a list of all genotypes, \code{mom},
#' mom's haplotypes, and \code{which} a dataframe of relatedness type and which
#' haplotypes an individual received.
sibFamily <- function(n, nsites, popsize=1000, selfing=0.5, fullsib=0.3, nfullsib=2,
                      mom_inbred=FALSE, ehet=0.8, ehom=0.1) {
  # create the population of grandparents of the current pop, which we use to
  # create gametes for the currrent pop.
  grandpop <- createDiploidSample(popsize, nsites)

  # shuffle parents
  p1 <- sample(1:popsize, popsize)
  p2 <- sample(1:popsize, popsize)

  # create gametes from this permutation
  pop <- mapply(function(m, f) {
                ma <- sampleAllele(grandpop[, m])
                fa <- sampleAllele(grandpop[, f])
                cbind(ma, fa)
  }, p1, p2, SIMPLIFY=FALSE)

  # make sibling family
  # events are the type of individual this is (selfed, fullsib, halfsib)
  ind_type <- sample(1:3, n, prob=c(selfing, fullsib, 1-(selfing+fullsib)),
                     replace=TRUE)

  # generate gametes according to these individual types
  # This is done in parts, first with halfsibs (these are most restrictive)
  hs_fathers <- sample(1:popsize, sum(ind_type == 3), replace=FALSE)
  fs_potential_fathers <- sample(setdiff(1:popsize, hs_fathers), 2, replace=FALSE)
  fs_fathers <- sample(fs_potential_fathers, sum(ind_type == 2), replace=TRUE)
  mom_i <- sample(setdiff(1:popsize, c(hs_fathers, fs_fathers)), 1)

  if (mom_inbred) {
    # make her gametes the same (choosing one randomly)
    i <- sample(1:2, 1)
    pop[[mom_i]] <- pop[[mom_i]][, c(i, i)]
  }

  # now mate individuals' haplotypes, taking a random haplotype from fathers
  hs <- mapply(function(f) sampleHaplotypes(pop[[mom_i]], pop[[f]]), hs_fathers,
               SIMPLIFY=FALSE)
  fs <- mapply(function(f) sampleHaplotypes(pop[[mom_i]], pop[[f]]), fs_fathers,
               SIMPLIFY=FALSE)
  # selfs:
  sf <- mapply(function(x) {
               sampleHaplotypes(pop[[mom_i]], pop[[mom_i]])
               }, seq_len(sum(ind_type == 1)), SIMPLIFY=FALSE)

  fs <- haplotypeSample2List(fs, "full sib", ehet, ehom)
  hs <- haplotypeSample2List(hs, "half sib", ehet, ehom)
  sf <- haplotypeSample2List(sf, "self", ehet, ehom)

  list(geno=cbind(hs$geno, fs$geno, sf$geno),
       mom=pop[[mom_i]],
       which=rbind(hs$which, fs$which, sf$which))
}


