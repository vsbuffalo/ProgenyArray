# parentage.R -- 
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

# not testing over this parameters space
NPARENT <- 100
NPROGENY <- 100
SELFING <- 0.5

error_message <- function(ex, act) sprintf("proportion correct below expectations: %0.3f expected, %0.3f actual", ex, act)

# test parentage, expecting a certain degree of correct calls for some set of parameters
test_parentage <- function(max_error, nloci, prop_parent_missing, prop_progeny_missing,
                           ehet, ehom, ehet_inf=ehet, ehom_inf=ehom) {
  set.seed(0)
  shutup <- function(x) suppressMessages(suppressWarnings(x))
  pa <- shutup(SimulatedProgenyArray(nparent=NPARENT, nprogeny=NPROGENY, nloci=nloci, selfing=SELFING,
                                     prop_parent_missing=prop_parent_missing, prop_progeny_missing=prop_progeny_missing,
                                     ehet=ehet, ehom=ehom))

  pa <- shutup(inferParents(pa, ehet_inf, ehom_inf))
  error <- propCorrectParents(pa)
  msg <- sprintf("parentage error rate exceeded with loci=%d, error=(%0.2f, %0.2f), missing=(%0.2f, %0.2f)",
                 nloci, ehet, ehom, prop_parent_missing, prop_progeny_missing)
  test_that(msg, {
            expect_true(error > max_error, error_message(max_error, error))})
}



## To regression test, we use simulated data and send bounds for errors
context("parentage tests, no error")
test_parentage(0.90, 50, 0, 0, 0, 0)
test_parentage(0.96, 100, 0, 0, 0, 0)

context("parentage tests, with error")
# 500 loci with error = (het=0.5, hom=0.1)
test_parentage(0.85, 500, 0, 0, 0.5, 0.1)

# 1000 loci with error = (het=0.5, hom=0.1)
test_parentage(0.98, 1000, 0, 0, 0.5, 0.1)

# 5000 loci with error = (het=0.5, hom=0.1) and missing = (0.01, 0.1)
test_parentage(0.95, 5000, 0.01, 0.1, 0.5, 0.1)
