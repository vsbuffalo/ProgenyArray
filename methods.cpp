#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector countGenotypes(IntegerVector x) {
	// This method is a specialized version of R's table that counts genotypes
	// encoded as 0, 1, 2 in a vector (and also returns NA) always of length 4,
	// always as numbers of 0, 1, 2, NA. This allows faster usage with apply, as
	// we don't need to convert to factor to get all genotype counts, even if
	// none are present.
	CharacterVector names = CharacterVector::create("0", "1", "2", "NA");
	IntegerVector genocounts(4);
	genocounts.fill(0);
	for (int i = 0; i < x.length(); i++) {
		if (!IntegerVector::is_na(x[i])) {
			genocounts[x[i]] += 1;
		} else {
			genocounts[3] += 1;
		}
	}
	genocounts.attr("names") = names;
	return genocounts;
}
