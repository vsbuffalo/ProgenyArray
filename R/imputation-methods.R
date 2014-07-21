# imputation-methods.R -- parent genotype imputation
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

parentGenotypeLikelihood <- function(halfsib_geno, freq, error_matrix) {
  # for single locus
  cp <- conditionalOffspringParentMatrix(freq)

  # outer loop: calculate the probability of that particular offspring
  #ccases <- complete.cases(halfsib_geno)
  #genos <- halfsib_geno[ccases, , drop=FALSE]
  genos <- halfsib_geno
  oprobs <- lapply(genos, function(g) {
                   # inner loop: calculate for a particular parent
                   sapply(1:3, function(i) {
                          # TODO the counts table should be a sufficient statistic, so this is unneccesary.
                          sum(error_matrix[, g+1L] * cp[i, ])
                   })
  })
  stopifnot(length(oprobs) > 0)
  do.call(rbind, oprobs)
}

parentGenotypesLikelihoods <- function(progeny_geno, freqs, ehet=0.8, ehom=0.1) {
  error_matrix <- genotypingErrorMatrix(ehet, ehom)
  mm <- cbind(freqs, progeny_geno)
  tmp <- apply(mm, 1, function(both) {
               f <- both[1]
               genos <- both[2:length(both)]
               hw <- c((1-f)^2, 2*f*(1-f), f^2)
               genos_cc <- na.omit(genos)
               if (!length(genos_cc))
                 m <- matrix(c(1, 1, 1), nrow=1)
               else
                 m <- parentGenotypeLikelihood(genos_cc, f, error_matrix)
               colSums(log(m)) + log(hw)
  })
  t(tmp)
}

which_max_NA <- function(x) {
  # which.max() that returns NA if all values areNA
  if (all(is.na(x)))
    return(NA)
  which.max(x)
}

imputeParentGenotypes <- function(x, freqs, ehet, ehom) {
#  i <- complete.cases(x)
#  if (length(i) == nrow(x)) stop("no complete cases")
#  if (length(i) < nrow(x))
#    warning(sprintf("%d loci not imputed due to completely missing cases in progeny", nrow(x) - length(i)))
#  x <- x[i, , drop=FALSE]
#  freqs <- freqs[i]
  ll <- parentGenotypesLikelihoods(x, freqs, ehet, ehom)
  ml_geno <- apply(ll, 1, which_max_NA) - 1L
  ml_ll <- apply(ll, 1, max)
  data.frame(geno=ml_geno, ll=ml_ll)
}

.getHalfsibFamilyGeno <- function(mother, mothers, fathers, progeny_geno, allowfs=FALSE) {
  #g <- progenyGenotypes(x)[, mothers(x) == mother]
  ## use the real mother: no point in adding uncertainty of parentage in this step
  ## we use ProgenyArray genotypes, as these have error, missing according to sim parameters
  g <- progeny_geno[, mothers == mother]
  # get the fathers of all kids from this mother
  dads <- fathers[mothers == mother]
  # ditch selfs
  not_selfed_or_full_sib <- dads != mother
  if (!allowfs) # ditch full sibs if not allowed
   not_selfed_or_full_sib <- not_selfed_or_full_sib & !duplicated(dads)
  g[, not_selfed_or_full_sib, drop=FALSE]
}

# TODO MOVE
setMethod("getHalfsibFamilyGeno", "ProgenyArray", function(x, mother) {
  .getHalfsibFamilyGeno(mother, mothers(x), fathers(x), progenyGenotypes(x))
})

setMethod("getHalfsibFamilyGeno", "SimulatedProgenyArray", function(x, mother) {
  .getHalfsibFamilyGeno(mother, x@mothers, x@fathers, x@progeny_geno)
})

