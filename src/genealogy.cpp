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
std::array<Matrix3d,3>
MendelianTransmissionMatrix() {
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

// Conditional offpsring given parent, P(G_o | G_p)
Matrix3d
ConditionalTransmissionMatrix(double freq) {
  Matrix3d tmprob;
  tmprob << 1 - freq,       freq,       0,
            (1 - freq)/2.,  1/2.,       freq/2.,
            0,              1 - freq,   freq;
  return tmprob;
}

// Genotype error matrix
Matrix3d
GenotypeError(double ehet, double ehom) {
  Matrix3d error;
  error << 1-ehom, ehom/2, ehom/2,
           ehet/2, 1-ehet, ehet/2,
           ehom/2, ehom/2, 1-ehom;
  return error;
}

// Product of Mendelian transmission matrix and genotype error
std::array<Matrix3d,3>
MendelianTransmissionWithErrorMatrix(double ehet, double ehom, bool use_log) {
  std::array<Matrix3d,3> tmprobs = MendelianTransmissionMatrix();
  Matrix3d error = GenotypeError(ehet, ehom);
  for (int i = 0; i < 3; i++) {
    if (use_log)
      tmprobs[i] = (tmprobs[i] * error).array().log().matrix();
    else
      tmprobs[i] = tmprobs[i] * error;
  }
  return tmprobs;
}

// Product of conditional offspring matrix and genotype error
Matrix3d
ConditionalTransmissionMatrixWithError(double freq, double ehet, double ehom, bool use_log) {
  Matrix3d tmprob = ConditionalTransmissionMatrix(freq);
  Matrix3d error = GenotypeError(ehet, ehom);
  if (use_log)
    tmprob = (tmprob * error).array().log().matrix();
  else
    tmprob = tmprob * error;
  return tmprob;
}


// [[Rcpp::export(".printMendelianMatrices")]]
void
printMendelianMatrices(NumericVector ehetv, NumericVector ehomv, LogicalVector log) {
  bool use_log = log[0];
  double ehet = ehetv[0], ehom = ehomv[0];
  Matrix3d tmprobs;
  //tmprobs = MendelianTransmissionWithErrorMatrix(ehet, ehom, use_log);
  tmprobs = ConditionalTransmissionMatrixWithError(0.5, ehet, ehom, use_log);
  //for (int i = 0; i < 3; i++)
  //  std::cout << tmprobs[i] << std::endl << std::endl;
  std::cout << tmprobs << std::endl << std::endl;
}

// [[Rcpp::export(".probQQ")]]
NumericVector
probQQ(IntegerVector progeny, IntegerVector parent_1, IntegerVector parent_2) {
//              std::array<Matrix3d,3> &transm) {
  // TODO size assertion
  std::array<Matrix3d,3> tmprobs;
  tmprobs = MendelianTransmissionWithErrorMatrix(0.8, 0.1, true);
  NumericVector out(1);
  double ll = 0;
  for (int i = 0; i < progeny.size(); i++) {
    ll += tmprobs[parent_1[i]](parent_2[i], progeny[i]);
  }
  out[0] = ll;
  return out;
}

// Probability of offspring given one parent and random allele draw from population
// [[Rcpp::export(".probQU")]]
NumericVector
probQU(NumericVector freqs, IntegerVector progeny, IntegerVector parent,
       double ehet, double ehom) {
  // TODO size assertion
  Matrix3d tmprob;
  NumericVector out(1);
  double ll = 0;
  for (int i = 0; i < progeny.size(); i++) {
    tmprob = ConditionalTransmissionMatrixWithError(freqs[i], ehet, ehom, true);
    ll += tmprob(parent[i], progeny[i]);
  }
  out[0] = ll;
  return out;
}

// Probability of seeing these genotypes from random allele draws in population
// [[Rcpp::export(".probUU")]]
NumericVector
probUU(NumericVector freqs, IntegerVector progeny, IntegerVector parent_1,
       IntegerVector parent_2, double ehet, double ehom) {
  Matrix3d tmprob;
  NumericVector out(1);
  double f, ll = 0;
  for (int i = 0; i < progeny.size(); i++) {
    f = freqs[i];
    std::array<double,3> hwprobs = {log(pow(1-f, 2.)),
                                    log(2*f*(1-f)),
                                    log(pow(f, 2.))};
    ll += hwprobs[progeny[i]];
    ll += hwprobs[parent_1[i]];
    ll += hwprobs[parent_2[i]];
  }
  out[0] = ll;
  return out;
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
