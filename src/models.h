#ifndef MODELS_H
#define MODELS_H


#include <Rcpp.h>
#include <RcppEigen.h>
#include <array>
#include <assert.h>
#include <vector>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::RowVector3d;

// Store the log probabilities under two relatedness models:
// two parents fully related (QQ) and unrelated (UU)
struct RelatednessProbs {
  double qq;
  double uu;
};

// Mendelian transmission probabilities are kept in same structure; make it a
// type
typedef std::array<Eigen::Matrix3d,3> MTMatrix;

MTMatrix MendelianTransmissionMatrix();
RowVector3d HWWithError(double freq, Matrix3d error, bool use_log);
Matrix3d ConditionalTransmissionMatrix(const double freq);
RowVector3d HWVector(double freq, bool use_log);
Matrix3d GenotypeError(const double ehet, const double ehom);

#endif /* MODELS_H */
