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

RowVector3d HWVector(double freq, bool use_log) {
  RowVector3d hwprobs(pow(1-freq, 2.), 2*freq*(1-freq), pow(freq, 2.));
  if (use_log)
    hwprobs = hwprobs.array().log().matrix();
  return hwprobs;
}

// TODO HERE const
RowVector3d HWWithError(double freq, Matrix3d error, bool use_log) {
  RowVector3d hwprobs(pow(1-freq, 2.), 2*freq*(1-freq), pow(freq, 2.));
  hwprobs = hwprobs * error;
  if (use_log)
    hwprobs = hwprobs.array().log().matrix();
  return hwprobs;
}

RelatednessProbs probQQUU(const MatrixXd &genofreqs,
    const MatrixXd &genofreqs_error,
    const MTMatrix &tmprobs, const IntegerVector progeny,
    const IntegerMatrix::Column parent_1, const IntegerMatrix::Column parent_2,
    const Matrix3d &errors, const LogicalVector use_locus) {
  double llqq = 0, lluu = 1;
  double m, f, p;
  RowVector3d hwprobs, log_hwprobs_error;
  for (int i = 0; i < progeny.length(); i++) { // for all loci
    if (!use_locus[i]) continue;
    hwprobs = genofreqs.row(i);
    log_hwprobs_error = genofreqs_error.row(i);
    // P(G | UU)
    lluu += log_hwprobs_error(progeny[i]) + log_hwprobs_error(parent_1[i]) + log_hwprobs_error(parent_2[i]);
    // P(G | QQ)
    // we need to sum over hidden genotypes of M, F, and O
    p = 0;
    for (int gm = 0; gm < 3; gm++) {
      for (int gf = 0; gf < 3; gf++) {
        for (int go = 0; go < 3; go++) {
          p += hwprobs[gf] * hwprobs[gm] * tmprobs[gm](gf, go) * errors(go, progeny[i]) * errors(gm, parent_1[i]) * errors(gf, parent_2[i]);
        }
      }
    }
    llqq += log(p);
  }
  RelatednessProbs out = {llqq, lluu};
  return out;
}

// allParentLikelihoods returns a list of each parent's log likelihoods under
// two models: both parents are related and neither are. This function also
// handles calculating a progeny's set of complete loci (e.g. not NA) that can
// be used for the parent likelihoods calculations
List parentageLikelihoods(const MatrixXd &genofreqs,
                          const MatrixXd &genofreqs_error,
                          const MTMatrix &tmprobs,
                          const IntegerVector progeny,
                          const IntegerMatrix parents,
                          const Matrix3d &errors) {
  int nparents = parents.ncol(), nloci = parents.nrow();
  NumericMatrix prob_qq(nparents, nparents), prob_uu(nparents, nparents);
  IntegerVector nloci_used(1);
  RelatednessProbs relprobs;

  // fill matrices with NA
  std::fill(prob_qq.begin(), prob_qq.end(), NumericMatrix::get_na());
  std::fill(prob_uu.begin(), prob_uu.end(), NumericMatrix::get_na());

  // A vector of non-NA progeny loci which is passed to other functions.
  // This prevents need for explicit subsetting and creating new vectors.
  // We also remove all fixed loci here.
  LogicalVector complete_loci(nloci);
  complete_loci = !is_na(progeny);
  nloci_used[0] = sum(complete_loci);

  for (int p1 = 0; p1 < nparents; p1++) {
    IntegerMatrix::Column par_1 = parents(_, p1);
    for (int p2 = p1; p2 < nparents; p2++) {
      IntegerMatrix::Column par_2 = parents(_, p2);
      relprobs = probQQUU(genofreqs, genofreqs_error, tmprobs,
                          progeny, par_1, par_2, errors, complete_loci);
      prob_qq(p1, p2) = relprobs.qq;
      prob_uu(p1, p2) = relprobs.uu;
    }
  }
  return List::create(_["prob_qq"]=prob_qq, _["prob_uu"]=prob_uu,
                      _["nloci"]=nloci_used);
}

// Calculate the HW genotype frequencies for all loci
MatrixXd HWLociMatrix(const NumericVector freqs, const bool use_log) {
  MatrixXd genofreqs(freqs.length(), 3);
  for (int i = 0; i < freqs.length(); i++) {
    genofreqs.row(i) = HWVector(freqs[i], use_log);
  }
  return genofreqs;
}

// Calculate the HW genotype frequencies for all loci
MatrixXd HWLociMatrixWithError(const NumericVector freqs,
    Matrix3d error, bool use_log) {
  MatrixXd genofreqs(freqs.length(), 3);
  for (int i = 0; i < freqs.length(); i++) {
    genofreqs.row(i) = HWVector(freqs[i], false) * error;
  }
  if (use_log)
    genofreqs = genofreqs.array().log().matrix();
  std::cout << genofreqs << std::endl;
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
  MatrixXd genofreqs, genofreqs_error;
  int nprogeny = progeny.ncol();
  List calls;
  bool verbose = verbosev[0];

  Matrix3d errors = GenotypeError(ehet, ehom);
  tmprobs = MendelianTransmissionMatrix();
  genofreqs = HWLociMatrix(freqs, false);
  genofreqs_error = HWLociMatrixWithError(freqs, errors, true);
  // TODOFIX
  for (int p = 0; p < nprogeny; p++) {
    checkUserInterrupt();
    IntegerMatrix::Column kid = progeny(_, p);
    if (verbose)
      std::cout << "\t" << p << "/" << nprogeny << " progeny completed\r" << std::flush;
    calls.push_back(parentageLikelihoods(genofreqs, genofreqs_error,
                                         tmprobs, kid, parents, errors));
  }
  if (verbose)
    std::cout << std::endl;
  return calls;
}


