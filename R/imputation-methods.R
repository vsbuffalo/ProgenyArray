# imputation-methods.R -- parent genotype imputation
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

imputeParentGenotype <- function(halfsib_geno, freq, error_matrix) {
  # for single locus
  cp <- conditionalOffspringParentMatrix(freq)

  # outer loop: calculate the probability of that particular offspring
  oprobs <- lapply(na.omit(halfsib_geno), function(g) {
                   # inner loop: calculate for a particular parent
                   sapply(1:3, function(i) {
                          sum(error_matrix[, g+1L] * cp[i, ])
                   })
  })
  do.call(rbind, oprobs)
}

imputeParentGenotypes <- function(progeny_geno, freqs, ehet=0.8, ehom=0.1) {
  error_matrix <- genotypingErrorMatrix(ehet, ehom)
  tmp <- apply(cbind(freqs, progeny_geno), 1, function(both) {
        f <- both[1]
        genos <- both[2:length(both)]
        hw <- c((1-f)^2, 2*f*(1-f), f^2)
        m <- imputeParentGenotype(genos, f, error_matrix)
        colSums(log(m)) + log(hw)
  })
  t(tmp)
}

getHalfsibFamilyGeno <- function(x, mothers, fathers, mother, allowfs=FALSE) {
  #g <- progenyGenotypes(x)[, mothers(x) == mother]
  g <- progenyGenotypes(x)[, mothers == mother]
  #dads <- fathers(x)[mothers(x) == mother]
  dads <- fathers[mothers == mother]
  not_selfed_or_full_sib <- dads != mother
  if (!allowfs)
   not_selfed_or_full_sib <- not_selfed_or_full_sib & !duplicated(dads)
  g[, not_selfed_or_full_sib, drop=FALSE]
}

# TODO genotype confusion matrix

