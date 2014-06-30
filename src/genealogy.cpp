#include <Rcpp.h>
#include <RcppEigen.h>
#include <array>
#include <assert.h>
#include <vector>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

// Note: frequencies vectors like `freq` always reference alternate allele
// frequency. Genotypes like 0, 1, 2, are ref/ref ref/alt, alt/alt.

using namespace Rcpp;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
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

// Product of Mendelian transmission matrix and genotype error
MTMatrix
MendelianTransmissionWithErrorMatrix(const double ehet, const double ehom,
                                     const bool use_log) {
  MTMatrix tmprobs = MendelianTransmissionMatrix();
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
ConditionalTransmissionMatrixWithError(const double freq, const double ehet,
                                       const double ehom, const bool use_log) {
  Matrix3d tmprob = ConditionalTransmissionMatrix(freq);
  Matrix3d error = GenotypeError(ehet, ehom);
  if (use_log)
    tmprob = (tmprob * error).array().log().matrix();
  else
    tmprob = tmprob * error;
  return tmprob;
}

// TODO HERE const
RowVector3d HWWithError(double freq, double ehet, double ehom, bool use_log) {
  RowVector3d hwprobs(pow(1-freq, 2.), 2*freq*(1-freq), pow(freq, 2.));
  hwprobs = hwprobs * GenotypeError(ehet, ehom);
  if (use_log)
    hwprobs = hwprobs.array().log().matrix();
  return hwprobs;
}

//// [[Rcpp::export(".printMendelianMatrices")]]
//void
//printMendelianMatrices(NumericVector ehetv, NumericVector ehomv, LogicalVector log) {
//  bool use_log = log[0];
//  double ehet = ehetv[0], ehom = ehomv[0];
//  Matrix3d tmprobs;
//  //tmprobs = MendelianTransmissionWithErrorMatrix(ehet, ehom, use_log);
//  tmprobs = ConditionalTransmissionMatrixWithError(0.5, ehet, ehom, use_log);
//  //for (int i = 0; i < 3; i++)
//  //  std::cout << tmprobs[i] << std::endl << std::endl;
//  std::cout << tmprobs << std::endl << std::endl;
//}

//double
//probQQ(const NumericVector freqs, const IntegerVector progeny,
//    const IntegerMatrix::Column parent_1, const IntegerMatrix::Column parent_2,
//    const double ehet, const double ehom, const std::array<Matrix3d,3> tmprobs,
//    const LogicalVector use_locus) {
//  // TODO size assertion
//  double ll = 0;
//  Vector3d hwprobs;
//  for (int i = 0; i < progeny.size(); i++) {
//    if (!use_locus[i]) continue;
//    hwprobs = HWWithError(freqs[i], ehet, ehom, true);
//    //std::cout << tmprobs[parent_1[i]] << std::endl << std::endl;
//    //if (progeny[i] > 2)
//    //  throw std::invalid_argument(""); // TODO
//    ll += tmprobs[parent_1[i]](parent_2[i], progeny[i]);
//    ll += hwprobs[parent_1[i]] + hwprobs[parent_2[i]];
//  }
//  return ll;
//}

// Probability of offspring given one parent and random allele draw from population
//double
//probQU(NumericVector freqs, IntegerVector progeny, IntegerMatrix::Column parent,
//       double ehet, double ehom, LogicalVector use_locus) {
//  // TODO size assertion
//  Matrix3d tmprob;
//  Vector3d hwprobs;
//  double ll = 0;
//  for (int i = 0; i < progeny.size(); i++) {
//    if (!use_locus[i]) continue;
//    hwprobs = HWWithError(freqs[i], ehet, ehom, true);
//    tmprob = ConditionalTransmissionMatrixWithError(freqs[i], ehet, ehom, true);
//    // TODO FIXME this probability calc needs to incorporate both parents
//    ll += tmprob(parent[i], progeny[i]);
//  }
//  return ll;
//}

// Probability of seeing these genotypes from random allele draws in population
//double
//probUU(const NumericVector freqs, const IntegerVector progeny,
//    const IntegerMatrix::Column parent_1, const IntegerMatrix::Column parent_2,
//    const double ehet, const double ehom, const LogicalVector use_locus) {
//  double ll = 0;
//  Vector3d hwprobs;
//  for (int i = 0; i < progeny.size(); i++) {
//    if (!use_locus[i]) continue;
//    hwprobs = HWWithError(freqs[i], ehet, ehom, true);
//    ll += hwprobs[progeny[i]];
//    ll += hwprobs[parent_1[i]];
//    ll += hwprobs[parent_2[i]];
//  }
//  return ll;
//}
//
// Calculate probability of both QQ and UU simultaneously (decreases number of
// iterations). This assumes observed paternal and maternal genotypes are correct
// in transmission, but spreads some probability around in genotype frequencies.
RelatednessProbs
probQQUU(const NumericVector freqs, const MatrixXd &genofreqs,
         const MTMatrix &tmprobs, const IntegerVector progeny,
         const IntegerMatrix::Column parent_1, const IntegerMatrix::Column parent_2,
         const double ehet, const double ehom, const LogicalVector use_locus) {
  double llqq = 0, lluu = 0;
  RowVector3d hwprobs;
  for (int i = 0; i < progeny.length(); i++) {
    if (!use_locus[i]) continue;
    hwprobs = genofreqs.row(i);
    // P(G | UU)
    lluu += hwprobs(progeny[i]);
    lluu += hwprobs(parent_1[i]);
    lluu += hwprobs(parent_2[i]);
    // P(G | QQ)
    llqq += tmprobs[parent_1[i]](parent_2[i], progeny[i]);
    llqq += hwprobs(parent_1[i]) + hwprobs(parent_2[i]);
  }
  RelatednessProbs out = {llqq, lluu};
  return out;
}


// allParentLikelihoods returns a list of each parent's log likelihoods under
// two models: both parents are related and neither are. This function also
// handles calculating a progeny's set of complete loci (e.g. not NA) that can
// be used for the parent likelihoods calculations
List parentageLikelihoods(const NumericVector freqs,
                          const MatrixXd &genofreqs,
                          const MTMatrix &tmprobs,
                          const IntegerVector progeny,
                          const IntegerMatrix parents,
                          const double ehet,
                          const double ehom) {
  int nparents = parents.ncol(), nloci = parents.nrow();
  NumericMatrix prob_qq(nparents, nparents), prob_uu(nparents, nparents);
  IntegerVector nloci_used(1);
  RelatednessProbs relprobs;

  // fill matrices with NA
  std::fill(prob_qq.begin(), prob_qq.end(), NumericMatrix::get_na());
  std::fill(prob_uu.begin(), prob_uu.end(), NumericMatrix::get_na());

  // A vector of non-NA progeny loci which is passed to other functions.
  // This prevents need for explicit subsetting and creating new vectors.
  LogicalVector complete_loci(nloci);
  complete_loci = !is_na(progeny);
  nloci_used[0] = sum(complete_loci);

  for (int p1 = 0; p1 < nparents; p1++) {
    IntegerMatrix::Column par_1 = parents(_, p1);
    for (int p2 = p1; p2 < nparents; p2++) {
      IntegerMatrix::Column par_2 = parents(_, p2);
      relprobs = probQQUU(freqs, genofreqs, tmprobs, progeny, par_1, par_2,
                          ehet, ehom, complete_loci);
      prob_qq(p1, p2) = relprobs.qq;
      prob_uu(p1, p2) = relprobs.uu;
    }
  }
  return List::create(_["prob_qq"]=prob_qq, _["prob_uu"]=prob_uu,
                      _["nloci"]=nloci_used);
}

// Calculate the HW genotype frequencies with error for all loci, with error
MatrixXd HWWithErrorMatrix(const NumericVector freqs, double ehet, double ehom,
                           const bool use_log) {
  MatrixXd genofreqs(freqs.length(), 3);
  for (int i = 0; i < freqs.length(); i++) {
    genofreqs.row(i) = HWWithError(freqs[i], ehet, ehom, use_log);
  }
  return genofreqs;
}

// Calculate the likelihoods of each parent under different models of
// relatedness. This requires that the progeny and parents matrices have
// already been subset so that all parental loci are complete (no NAs). NAs in
// the kid matrix are allowed, as allParentLikelihoods will only operate on the
// subset of loci the kid has complete.
// [[Rcpp::export(".inferParents")]]
List inferParents(IntegerMatrix progeny, IntegerMatrix parents, NumericVector
    freqs, NumericVector ehetv, NumericVector ehomv, LogicalVector verbosev) {

  double ehet = ehetv[0], ehom = ehomv[0];
  MTMatrix tmprobs;
  MatrixXd genofreqs;
  int nprogeny = progeny.ncol();
  List calls;
  bool verbose = verbosev[0];

  tmprobs = MendelianTransmissionWithErrorMatrix(ehet, ehom, true);
  genofreqs = HWWithErrorMatrix(freqs, ehet, ehom, true);

  for (int p = 0; p < nprogeny; p++) {
    checkUserInterrupt();
    IntegerMatrix::Column kid = progeny(_, p);
    if (verbose)
      std::cout << "\t" << p << "/" << nprogeny << " progeny completed\r" << std::flush;
    calls.push_back(parentageLikelihoods(freqs, genofreqs, tmprobs, kid, parents, ehet, ehom));
  }
  if (verbose)
    std::cout << std::endl;
  return calls;
}


