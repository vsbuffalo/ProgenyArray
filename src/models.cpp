#include <Rcpp.h>
#include <RcppEigen.h>
#include <array>
#include <assert.h>
#include <vector>
#include <iostream>
#include "models.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::RowVector3d;
using Eigen::Vector2d;
using Eigen::Array2d;
using Eigen::Vector3d;



// Mendelian transmission matrix, P(G_o | G_m, G_f)
// dimensions: mother, father, offspring
MTMatrix
MendelianTransmissionMatrix() {
  double transprobs[3][9] = {
    {1., 1/2., 0, 0, 1/2., 1., 0, 0, 0},
    {1/2., 1/4., 0, 1/2., 1/2., 1/2., 0, 1/4., 1/2.},
    {0, 0, 0, 1., 1/2., 0, 0, 1/2., 1.}};
  MTMatrix tmprobs;
  for (int i = 0; i < 3; i++) {
    tmprobs[i] = Map<Matrix3d>(transprobs[i]);
  }
  return tmprobs;
}

// Conditional offspring given parent, P(G_o | G_p)
Matrix3d
ConditionalTransmissionMatrix(const double freq) {
  Matrix3d tmprob;
  tmprob << 1 - freq,       freq,       0,
            (1 - freq)/2.,  1/2.,       freq/2.,
            0,              1 - freq,   freq;
  return tmprob;
}

// Genotype error matrix
Matrix3d
GenotypeError(const double ehet, const double ehom) {
  Matrix3d error;
  error << 1-ehom, ehom/2, ehom/2,
           ehet/2, 1-ehet, ehet/2,
           ehom/2, ehom/2, 1-ehom;
  return error;
}

RowVector3d HWVector(double freq, bool use_log) {
  RowVector3d hwprobs(pow(1-freq, 2.), 2*freq*(1-freq), pow(freq, 2.));
  if (use_log)
    hwprobs = hwprobs.array().log().matrix();
  return hwprobs;
}


