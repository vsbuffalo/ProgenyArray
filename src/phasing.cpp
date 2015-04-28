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
using Eigen::VectorXd;
using Eigen::Array2d;
using Eigen::Vector3d;

typedef Eigen::Matrix<double, 3, 2> ATMatrix;

ATMatrix AlleleTransmissionMatrix() {
  Eigen::Matrix<double, 3, 2> m;
  m << 1, 0, 1/2., 1/2., 0, 1;
  return m;
};

Vector2d pa_gp(const int gp, const double freq, Matrix3d error,
 RowVector3d hw, ATMatrix pag) {
  // P(a | g'), returns 2d row vector of allele probs given observed genotype
  // of parent.
  Vector2d probs(0, 0);
  double denom=0;
  int a;
  for (int i=0; i < 3; i++) denom += error(i, gp) * hw(i);
  for (int i=0; i < 3; i++) {
    for (a=0; a < 2; a++) {
      probs(a) += hw(i) * error(i, gp) * pag(i, a);
    }
  }
  probs = (probs.array() / denom).matrix();
  return(probs);
}

// [[Rcpp::export]]
Eigen::Vector2d pa_gp_ext(int gp, double freq, double ehet, double ehom) {
  // exported version (for debugging)
  RowVector3d hw = HWVector(freq, false);
  ATMatrix pag = AlleleTransmissionMatrix();
  Matrix3d error = GenotypeError(ehet, ehom);
  return pa_gp(gp, freq, error, hw, pag);
}

Eigen::Matrix<double, 2, 3> pg_theta(const int gp, const double freq,
    Matrix3d error, RowVector3d hw, ATMatrix pag) {
  // P(g | theta, g_f')
  Eigen::Matrix<double, 2, 3> m;
  Vector2d pagp = pa_gp(gp, freq, error, hw, pag);
  m << pagp(0), pagp(1),       0,
       0      , pagp(0), pagp(1);
  return(m);
}

// [[Rcpp::export]]
NumericMatrix pg_theta_ext(int gp, double freq, double ehet, double ehom) {
  // exported version (for debugging)
  RowVector3d hw = HWVector(freq, false);
  ATMatrix pag = AlleleTransmissionMatrix();
  Matrix3d error = GenotypeError(ehet, ehom);
  Eigen::Matrix<double, 2, 3> pgt;
  pgt = pg_theta(gp, freq, error, hw, pag);
  NumericMatrix coef(wrap(pgt));
  return coef;
}

//NumericVector allele_ll(IntegerVector progeny_genos, double ehet, double ehom,
//  // deprecated
//    double freq, IntegerVector father_genos) {
//  Vector2d ll(0, 0), p(0, 0);
//  RowVector3d hw = HWVector(freq, false);
//  ATMatrix pag = AlleleTransmissionMatrix();
//  Matrix3d error = GenotypeError(ehet, ehom);
//  for (int j=0; j < progeny_genos.length(); j++) {
//    if (IntegerVector::is_na(progeny_genos(j)) || IntegerVector::is_na(father_genos(j)))
//      continue;
//    p << 0, 0;
//    for (int i=0; i < 3; i++) {
//      // TODO we could refactor to use lin algebra
//      p += error(i, progeny_genos(j)) * pg_theta(father_genos(j), freq, error, hw, pag).col(i);
//    }
//    ll += p.array().log().matrix();
//  }
//  NumericVector llout(wrap(ll));
//  return llout;
//}

// [[Rcpp::export]]
NumericMatrix max_ll_matrix(NumericMatrix x, NumericMatrix y, IntegerVector which_max) {
  NumericMatrix out(x.nrow(), x.ncol());
  for (int i=0; i < x.nrow(); i++) {
    if (which_max[i] == 0)
      out.row(i) = x.row(i);
    else if (which_max[i] == 1)
      out.row(i) = y.row(i);
    else
      std::cout << "ERROR: which_max values must be in (0, 1)" << std::endl;
  }
  return out;
}

// [[Rcpp::export]]
List geno_ll(IntegerMatrix progeny_genos, IntegerMatrix father_genos,
    IntegerVector which_father, NumericVector freqs, double ehet, double ehom) {
  // P(geno | theta, father_geno)
  // for EM algorithm, M-step
  //VectorXd ll(progeny_genos.nrow(), 2);
  NumericMatrix ll0(progeny_genos.nrow(), progeny_genos.ncol());
  NumericMatrix ll1(progeny_genos.nrow(), progeny_genos.ncol());
  List lls;
  RowVector3d hw;
  Vector2d p;
  Eigen::Matrix<double, 2, 3> tmp;
  ATMatrix pag = AlleleTransmissionMatrix(); // TODO move this out
  Matrix3d error = GenotypeError(ehet, ehom); // TODO move this out
  double ll, pi;
  int this_father;

  if (progeny_genos.nrow() != father_genos.nrow()) {
    std::cout << "ERROR: unequal number of progeny and father loci!" << std::endl;
    return R_NilValue;
  }

  for (int i=0; i < progeny_genos.ncol(); i++) {
    // for each individual, calculate r_k * p(x_i | \theta_k)
    this_father = which_father[i];
    if (this_father >= father_genos.nrow()) {
      std::cout << "ERROR: this father's index (" << this_father << ") outside of bounds of father genotype matrix (" << father_genos.nrow() << ")" << std::endl;
      return R_NilValue;
    }
    for (int l=0; l < progeny_genos.nrow(); l++) {
      if (IntegerVector::is_na(progeny_genos(l, i)) ||
          IntegerVector::is_na(father_genos(l, this_father))) {
        ll0(l, i) = NA_REAL;
        ll1(l, i) = NA_REAL;
        continue;
      }
      hw = HWVector(freqs[l], false);
      tmp = pg_theta(father_genos(l, this_father), freqs[l], error, hw, pag) * error;
      p = tmp.col(progeny_genos(l, i));
      ll0(l, i) = p[0];
      ll1(l, i) = p[1];
    }
  }
  lls.push_back(ll0);
  lls.push_back(ll1);
  return lls;
}
  // TODO
  // freqs x
  // NA x
  // father indices x
  //
  // TODO sanity check:
  // no NA in father matrix
  // no complete loci with NA


