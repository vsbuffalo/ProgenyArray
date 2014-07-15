# ProgenyArray

[![Build Status](https://travis-ci.org/vsbuffalo/ProgenyArray.svg?branch=master)](https://travis-ci.org/vsbuffalo/ProgenyArray)

ProgenyArray includes classes and methods for working with genotyping data from
a progeny array.

This package implements some simple classes for storing data from a progeny
array (half-sib) experiment, as well as some statistical genetics methods to do
tasks like infer parentage and impute incorrect or missing genotype. These
methods were designed to be robust to very noisy data; the project this is
being applied to uses very low-coverage GBS data with about 40% data, 80%
heterozygous error rate, and a 10% homozygous error rate.


