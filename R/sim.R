## sim.R -- simulation of progeny array data


createDiploidSample  <- function(n, nsites) {
	# this is a bit of an approximation: create SFS, and then generate site
	# frequencies from that.

	sf <- 1/seq_len(n) # sfs for sample of n
	site_probs <- sf/sum(sf) # probability of getting a certain number of sites
	site_freqs <- seq_len(n)/n

	sample_sf <- sample(site_freqs, nsites, replace=TRUE, prob=site_probs)
	#samples <- t(replicate(2*n, { Vectorize(function(prob) rbinom(1, 1, prob), "prob")(sample_sf) } ))
	#split(split(samples, row(samples)), rep(seq_len(n), each=2))
  out <- lapply(seq_len(n), function(x) rbinom(nsites, 1, sample_sf)+rbinom(nsites, 1, sample_sf))
  do.call(cbind, out)
}

mate <- function(x, y) {
  # mate two genotype vectors
  # meosis
  transmit <- list(c(0), c(0, 1), c(1, 1))
  x_gametes <- sapply(transmit[x+1L], sample, size=1)
  y_gametes <- sapply(transmit[y+1L], sample, size=1)
  return(x_gametes + y_gametes)
}

createProgenyArrayWithSelfing <- function(x, ehet=0, ehom=0, nprogeny=50) {
  # create progeny array for a single mother from a population x
  mom_i <- sample(seq_len(ncol(x)), 1)
  dad_i <- sample(seq_len(ncol(x)), nprogeny, replace=TRUE)
  genotype_error <- genotypingErrorMatrix(ehet, ehom)
  progeny <- do.call(cbind, lapply(dad_i, function(di) {
                                   o <- mate(x[, mom_i], x[, di])
                                   sapply(o, function(y) sample(c(0L, 1L, 2L), 1, prob=genotype_error[y+1L, ]))
                                   }))
  list(progeny=progeny, mom_i=mom_i, dad_i=dad_i)
}

#createProgenyArrayWithSelfing <- function(x, nprogeny=50) {
#  # create progeny array for a single mother from a population x
#  mom_i <- sample(seq_len(ncol(x)), 1)
#  dad_i <- sample(seq_len(ncol(x)), nprogeny, replace=TRUE)
#  progeny <- do.call(cbind, lapply(dad_i, function(di) mate(x[, mom_i], x[, di])))
#  list(progeny=progeny, mom_i=mom_i, dad_i=dad_i)
#}
#

## testing
tm <- transmissionMatrix()
x <- createDiploidSample(100, 1000)

# with no errors
prog <- createProgenyArrayWithSelfing(x)
correct <- logical(length(prog$dad_i))
results <- vector("list", length(prog$dad_i))
mother <- x[, prog$mom_i]
for (i in seq_len(ncol(prog$progeny))) {
  res <- inferFather(prog$progeny[, i], mother, x, 0.1, 0.1)
  results[[i]] <- res
  correct[i] <- res$father == prog$dad_i[i]
}

# with errors
prog <- createProgenyArrayWithSelfing(x, 0.6, 0.1)
correct <- logical(length(prog$dad_i))
results <- vector("list", length(prog$dad_i))
mother <- x[, prog$mom_i]
for (i in seq_len(ncol(prog$progeny))) {
  res <- inferFather(prog$progeny[, i], mother, x, 0.6, 0.1)
  results[[i]] <- res
  correct[i] <- res$father == prog$dad_i[i]
}

