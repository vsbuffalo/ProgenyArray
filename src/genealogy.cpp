#include <Rcpp.h>
#include <RcppEigen.h>
#include <array>
#include <vector>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::Matrix3d;


// Mendelian transmission matrix
std::array<Matrix3d,3> MendelianMatrices() {
  double transprobs[3][9] = {
    {1., 1/2., 0, 0, 1/2., 1., 0, 0, 0},
    {1/2., 1/4., 0, 1/2., 1/2., 1/2., 0, 1/4., 1/2.},
    {0, 0, 0, 1., 1/2., 0, 0, 1/2., 1.}};
  std::array<Matrix3d,3> tmprobs;
  for (int i = 0; i < 3; i++) {
    tmprobs[i] = Map<Matrix3d>(transprobs[i]);
  }
  return tmprobs;
}

// Genotype error matrix
Matrix3d GenotypeError(double ehet, double ehom) {
  Matrix3d error;
  error << 1-ehom, ehom/2, ehom/2,
           ehet/2, 1-ehom, ehet/2,
           ehom/2, ehom/2, 1-ehom;
  return error;
}

// Product of Mendelian transmission matrix and genotype error
std::array<Matrix3d,3>
MendelianMatricesWithError(double ehet, double ehom, bool use_log) {
  std::array<Matrix3d,3> tmprobs = MendelianMatrices();
  Matrix3d error = GenotypeError(ehet, ehom);
  for (int i = 0; i < 3; i++) {
    if (use_log)
      tmprobs[i] = (tmprobs[i] * error).array().log().matrix();
    else
      tmprobs[i] = tmprobs[i] * error;
  }
  return tmprobs;
}

// [[Rcpp::export(".printMendelianMatrices")]]
void printMendelianMatrices(NumericVector ehetv, NumericVector ehomv, LogicalVector log) {
  bool use_log = log[0];
  double ehet = ehetv[0], ehom = ehomv[0];
  std::array<Matrix3d,3> mendel_matrix;
  mendel_matrix = MendelianMatricesWithError(ehet, ehom, use_log);
  for (int i = 0; i < 3; i++)
    std::cout << mendel_matrix[i] << std::endl << std::endl;
}

// [[Rcpp::export(".probQQ")]]
double probQQ(IntegerVector &progeny, IntegerVector &parent_1, IntegerVector &parent_2,
              std::array<Matrix3d,3> &transm) {
  // TODO size assertion
  double ll = 0;
  for (i = 0; i < progeny.size(); i++) {
    ll += transm[parent_1[i]][parent_2[i]][progeny[i]];
  }
  return ll;
}


//// [[Rcpp::export(".genealogyLikelihood")]]

//List maxLikelihoodGenealogy(IntegerMatrix parents, NumericVector freqs, NumericVector ehetv, NumericVector ehomv) {
//  // Find the maximum likelihood genealogy for a trio (allowing for selfing
//  // too). For P parents, this evaluates:
//  //
//  // (1) the probability of each of the P(P - 1)/2 possible two parents being
//  // parents
//  // (2) the probability of each of the 20 parents and a random parent
//  // from the population (there are P of these)
//  // (3) the probabiity of each parent having selfed (P of these)

//  // (4) the probablity of two random parents from the population
//  //
//  int nloci = parents.nrow(), nparents = parents.ncol();
//  double ehet = ehetv[0], ehom=ehom[0];
//
//  // Allocate a matrix for all parent comparisons for (1)
//  NumericMatrix prob_qq(nparents, nparents);
//  // Allocate vectors for each of the one-related parent, selfing, and random parent probabilties
//  NumericVector prob_qu(nparents);
//  NumericVector prob_q(nparents);
//  NumericVector prob_uu(nparents);
//
//  for (int p1 = 0; p < nparents; p++) {
//    for (int p2 = 0; p < nparents; p++) {
//
//    }
//  }
//
//}
